"""
UniProtDB — query interface for the local UniProt reference database.

Provides bidirectional lookup:
  sequence  →  similarity hits  →  domain/family annotations
  domain/family  →  representative sequences  →  cluster members

Quick start
-----------
    from uniprotdb import UniProtDB

    db = UniProtDB(
        db_path="/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/uniprot.sqlite",
        mmseqs_db_path="/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/mmseqs_db/uniref90_mmseqs",
    )

    # Sequence similarity search
    hits = db.search_by_sequence(fasta_str)

    # Fetch a specific protein
    info = db.get_by_accession("P04637")

    # Domain-based lookup
    proteins = db.search_by_domain(pfam_id="PF00069")

    # Start the MMseqs2 server before heavy usage:
    #   bash mmseqs_server.sh
"""

import io
import logging
import os
import shutil
import socket
import sqlite3
import subprocess
import time
import uuid
from contextlib import contextmanager
from datetime import date
from typing import Dict, List, Optional

import pandas as pd
import requests

log = logging.getLogger(__name__)

# UniProt REST API
_UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
_RATE_LIMIT_SLEEP  = 0.34          # ≤3 req/s
_SERVER_PORT       = 8080
_SERVER_URL        = f"http://localhost:{_SERVER_PORT}"


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _open_db(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path, check_same_thread=False)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA cache_size=-2000000")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.row_factory = sqlite3.Row
    return conn


def _query_len(fasta_str: str) -> int:
    """Return total residue count from a FASTA string (strips header lines)."""
    return sum(
        len(line)
        for line in fasta_str.splitlines()
        if line and not line.startswith(">")
    )


def _parse_m8(tsv: str) -> pd.DataFrame:
    """Parse MMseqs2 m8-format result TSV into a DataFrame."""
    cols = [
        "query", "target", "pident", "alnlen", "mismatch",
        "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore",
    ]
    if not tsv.strip():
        return pd.DataFrame(columns=cols)
    return pd.read_csv(io.StringIO(tsv), sep="\t", names=cols, comment="#")


# ---------------------------------------------------------------------------
# UniProtDB class
# ---------------------------------------------------------------------------

class UniProtDB:
    """
    Bidirectional query interface for the local UniProt reference database.

    Parameters
    ----------
    db_path : str
        Path to uniprot.sqlite
    mmseqs_db_path : str
        Path to the MMseqs2 database prefix (e.g. …/mmseqs_db/uniref90_mmseqs)
    """

    def __init__(self, db_path: str, mmseqs_db_path: str) -> None:
        self.db_path        = db_path
        self.mmseqs_db_path = mmseqs_db_path
        self._conn: Optional[sqlite3.Connection] = None

    # -----------------------------------------------------------------------
    # Internal
    # -----------------------------------------------------------------------

    @contextmanager
    def _get_conn(self):
        """Yield a SQLite connection; creates a persistent one on first call."""
        if self._conn is None:
            self._conn = _open_db(self.db_path)
        try:
            yield self._conn
        except Exception:
            self._conn.rollback()
            raise

    def close(self) -> None:
        if self._conn:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    # -----------------------------------------------------------------------
    # Server management
    # -----------------------------------------------------------------------

    def server_status(self) -> Dict:
        """
        Returns dict: {running, port, index_path, uptime}.

        Uses a TCP probe so no HTTP round-trip overhead when the server is down.
        uptime is None (the MMseqs2 server API does not expose it).
        """
        running = False
        try:
            with socket.create_connection(("localhost", _SERVER_PORT), timeout=2):
                running = True
        except (ConnectionRefusedError, OSError):
            pass

        return {
            "running":    running,
            "port":       _SERVER_PORT,
            "index_path": self.mmseqs_db_path,
            "uptime":     None,
        }

    # -----------------------------------------------------------------------
    # Sequence search
    # -----------------------------------------------------------------------

    def search_by_sequence(
        self,
        fasta_str: str,
        sensitivity: float = 7.5,
        max_hits: int = 100,
        min_identity: float = 0.3,
        min_coverage: float = 0.5,
    ) -> pd.DataFrame:
        """
        Search UniRef90 by sequence similarity.

        Returns
        -------
        DataFrame with columns:
          uniref90_id, rep_accession, evalue, identity, coverage,
          interpro_ids, pfam_ids, member_count, taxonomy_id
        """
        status = self.server_status()
        if status["running"]:
            tsv = self._search_via_server(fasta_str, sensitivity, max_hits, min_identity, min_coverage)
        else:
            log.info("MMseqs2 server not running — falling back to easy-search subprocess.")
            tsv = self._search_via_subprocess(fasta_str, sensitivity, max_hits, min_identity, min_coverage)

        return self._annotate_hits(tsv, fasta_str, min_identity, min_coverage)

    def _search_via_server(
        self,
        fasta_str: str,
        sensitivity: float,
        max_hits: int,
        min_identity: float,
        min_coverage: float,
    ) -> str:
        """POST to the running MMseqs2 REST server."""
        # MMseqs2 server accepts plain FASTA as request body.
        # Parameters are passed as query-string arguments.
        # NOTE: exact parameter names may vary with MMseqs2 version —
        #       verify against `mmseqs server --help` output.
        params = {
            "sensitivity": sensitivity,
            "max-seqs":    max_hits,
            "min-seq-id":  min_identity,
            "coverage":    min_coverage,
            "format-output": "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
        }
        resp = requests.post(
            f"{_SERVER_URL}/search",
            data=fasta_str.encode(),
            params=params,
            headers={"Content-Type": "text/plain"},
            timeout=600,
        )
        resp.raise_for_status()
        return resp.text

    def _search_via_subprocess(
        self,
        fasta_str: str,
        sensitivity: float,
        max_hits: int,
        min_identity: float,
        min_coverage: float,
    ) -> str:
        """Fallback: run mmseqs easy-search in a uniquely named temp directory."""
        run_id  = uuid.uuid4().hex
        tmp_dir = f"/tmp/mmseqs_{run_id}"
        os.makedirs(tmp_dir, exist_ok=True)

        query_file  = os.path.join(tmp_dir, "query.fasta")
        result_file = os.path.join(tmp_dir, "result.tsv")
        tmp_work    = os.path.join(tmp_dir, "work")
        os.makedirs(tmp_work)

        try:
            with open(query_file, "w") as fh:
                fh.write(fasta_str)

            cmd = [
                "mmseqs", "easy-search",
                query_file,
                self.mmseqs_db_path,
                result_file,
                tmp_work,
                "--format-output",
                "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
                "-s",          str(sensitivity),
                "--max-seqs",  str(max_hits),
                "--min-seq-id", str(min_identity),
                "-c",          str(min_coverage),
                "--threads",   "8",
            ]
            result = subprocess.run(
                cmd, check=True, capture_output=True, text=True
            )
            if result.stderr:
                log.debug("mmseqs stderr: %s", result.stderr[:500])

            if os.path.exists(result_file):
                with open(result_file) as fh:
                    return fh.read()
            return ""
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    def _annotate_hits(
        self,
        tsv: str,
        fasta_str: str,
        min_identity: float,
        min_coverage: float,
    ) -> pd.DataFrame:
        """Parse m8 results and join domain / cluster annotations from SQLite."""
        empty = pd.DataFrame(columns=[
            "uniref90_id", "rep_accession", "evalue", "identity",
            "coverage", "interpro_ids", "pfam_ids", "member_count", "taxonomy_id",
        ])

        df = _parse_m8(tsv)
        if df.empty:
            return empty

        q_len = _query_len(fasta_str) or 1

        df = df.rename(columns={"target": "uniref90_id", "pident": "identity_pct", "evalue": "evalue"})
        df["identity"] = df["identity_pct"] / 100.0
        # Query coverage: fraction of query aligned  (qstart/qend are 1-based)
        df["coverage"] = (df["qend"] - df["qstart"] + 1) / q_len

        df = df[df["identity"] >= min_identity]
        df = df[df["coverage"] >= min_coverage]
        if df.empty:
            return empty

        target_ids = df["uniref90_id"].tolist()

        with self._get_conn() as conn:
            # --- Sequence cluster metadata ---
            ph = ",".join("?" * len(target_ids))
            seq_rows = conn.execute(
                f"SELECT uniref90_id, rep_accession, member_count, taxonomy_id"
                f" FROM sequences WHERE uniref90_id IN ({ph})",
                target_ids,
            ).fetchall()
            seq_df = pd.DataFrame([dict(r) for r in seq_rows])

            if seq_df.empty:
                return empty

            df = df.merge(seq_df, on="uniref90_id", how="inner")

            # --- Domain annotations via rep_accession ---
            rep_accs = df["rep_accession"].dropna().unique().tolist()
            if rep_accs:
                ph2 = ",".join("?" * len(rep_accs))
                dom_rows = conn.execute(
                    f"""SELECT
                            uniprot_accession,
                            GROUP_CONCAT(DISTINCT CASE WHEN interpro_id IS NOT NULL
                                         THEN interpro_id END) AS interpro_ids,
                            GROUP_CONCAT(DISTINCT CASE WHEN db_name = 'Pfam'
                                         THEN db_accession END) AS pfam_ids
                        FROM domains
                        WHERE uniprot_accession IN ({ph2})
                        GROUP BY uniprot_accession""",
                    rep_accs,
                ).fetchall()
                dom_df = pd.DataFrame([dict(r) for r in dom_rows])
                dom_df = dom_df.rename(columns={"uniprot_accession": "rep_accession"})
                df = df.merge(dom_df, on="rep_accession", how="left")
            else:
                df["interpro_ids"] = None
                df["pfam_ids"]     = None

        out_cols = [
            "uniref90_id", "rep_accession", "evalue", "identity",
            "coverage", "interpro_ids", "pfam_ids", "member_count", "taxonomy_id",
        ]
        for c in out_cols:
            if c not in df.columns:
                df[c] = None

        return df[out_cols].sort_values("evalue").reset_index(drop=True)

    # -----------------------------------------------------------------------
    # Accession lookup
    # -----------------------------------------------------------------------

    def get_by_accession(self, uniprot_acc: str) -> Dict:
        """
        Full record for a UniProtKB accession.

        Returns
        -------
        dict with keys: sequence, domains, sprot_features, cluster_id,
                        is_reviewed, cluster_member_count
        """
        with self._get_conn() as conn:
            # Sequence + reviewed flag from cache
            cache_row = conn.execute(
                "SELECT sequence, reviewed FROM sequences_cache WHERE accession = ?",
                (uniprot_acc,),
            ).fetchone()

            if cache_row:
                sequence   = cache_row["sequence"]
                is_reviewed = bool(cache_row["reviewed"])
            else:
                sequence    = self.fetch_member_sequence(uniprot_acc)
                is_reviewed = False

            # Cluster membership
            id_row = conn.execute(
                "SELECT uniref90_id FROM idmapping WHERE uniprot_accession = ? LIMIT 1",
                (uniprot_acc,),
            ).fetchone()
            cluster_id = id_row["uniref90_id"] if id_row else None

            member_count = None
            if cluster_id:
                seq_row = conn.execute(
                    "SELECT member_count FROM sequences WHERE uniref90_id = ?",
                    (cluster_id,),
                ).fetchone()
                if seq_row:
                    member_count = seq_row["member_count"]

            # Domain annotations
            dom_rows = conn.execute(
                """SELECT interpro_id, interpro_name, db_name, db_accession,
                          start_pos, end_pos
                   FROM domains WHERE uniprot_accession = ?
                   ORDER BY start_pos""",
                (uniprot_acc,),
            ).fetchall()
            domains = [dict(r) for r in dom_rows]

            # Swiss-Prot features
            feat_rows = conn.execute(
                """SELECT feature_type, start_pos, end_pos, description
                   FROM sprot_features WHERE uniprot_accession = ?
                   ORDER BY start_pos""",
                (uniprot_acc,),
            ).fetchall()
            sprot_features = [dict(r) for r in feat_rows]

        return {
            "sequence":            sequence,
            "domains":             domains,
            "sprot_features":      sprot_features,
            "cluster_id":          cluster_id,
            "is_reviewed":         is_reviewed,
            "cluster_member_count": member_count,
        }

    # -----------------------------------------------------------------------
    # Domain-based search
    # -----------------------------------------------------------------------

    def search_by_domain(
        self,
        interpro_id: Optional[str] = None,
        pfam_id: Optional[str] = None,
        feature_type: Optional[str] = None,
        feature_description: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Find proteins annotated with a given domain or functional feature.

        At least one parameter must be provided.

        Returns
        -------
        DataFrame: uniprot_accession, uniref90_id, domain_info,
                   start_pos, end_pos, member_count
        """
        if not any([interpro_id, pfam_id, feature_type, feature_description]):
            raise ValueError("Provide at least one search criterion.")

        rows = []

        with self._get_conn() as conn:
            # --- InterPro / Pfam domains ---
            if interpro_id or pfam_id:
                sql  = """SELECT d.uniprot_accession, i.uniref90_id,
                                 d.db_name || ':' || COALESCE(d.db_accession,'') ||
                                 ' (' || COALESCE(d.interpro_id,'') || ') ' ||
                                 COALESCE(d.interpro_name,'') AS domain_info,
                                 d.start_pos, d.end_pos,
                                 s.member_count
                          FROM domains d
                          LEFT JOIN idmapping i ON i.uniprot_accession = d.uniprot_accession
                          LEFT JOIN sequences s ON s.uniref90_id = i.uniref90_id
                          WHERE """
                cond = []
                params: List = []
                if interpro_id:
                    cond.append("d.interpro_id = ?")
                    params.append(interpro_id)
                if pfam_id:
                    cond.append("(d.db_name = 'Pfam' AND d.db_accession = ?)")
                    params.append(pfam_id)
                sql += " OR ".join(cond)
                for r in conn.execute(sql, params).fetchall():
                    rows.append(dict(r))

            # --- Swiss-Prot features ---
            if feature_type or feature_description:
                sql  = """SELECT f.uniprot_accession, i.uniref90_id,
                                 f.feature_type || ': ' || COALESCE(f.description,'') AS domain_info,
                                 f.start_pos, f.end_pos,
                                 s.member_count
                          FROM sprot_features f
                          LEFT JOIN idmapping i ON i.uniprot_accession = f.uniprot_accession
                          LEFT JOIN sequences s ON s.uniref90_id = i.uniref90_id
                          WHERE """
                cond = []
                params = []
                if feature_type:
                    cond.append("f.feature_type = ?")
                    params.append(feature_type)
                if feature_description:
                    cond.append("f.description LIKE ?")
                    params.append(f"%{feature_description}%")
                sql += " AND ".join(cond)
                for r in conn.execute(sql, params).fetchall():
                    rows.append(dict(r))

        if not rows:
            return pd.DataFrame(columns=[
                "uniprot_accession", "uniref90_id", "domain_info",
                "start_pos", "end_pos", "member_count",
            ])
        return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)

    # -----------------------------------------------------------------------
    # Domain architecture
    # -----------------------------------------------------------------------

    def get_domain_architecture(self, uniprot_acc: str) -> List[Dict]:
        """
        Return all domain/feature annotations for a protein, sorted by start_pos.

        Returns
        -------
        list of dicts: {db_name, db_accession, interpro_id, name, start, end}
        """
        with self._get_conn() as conn:
            rows = conn.execute(
                """SELECT db_name, db_accession, interpro_id,
                          interpro_name AS name, start_pos AS start, end_pos AS end
                   FROM domains WHERE uniprot_accession = ?
                   ORDER BY start_pos""",
                (uniprot_acc,),
            ).fetchall()
        return [dict(r) for r in rows]

    # -----------------------------------------------------------------------
    # Cluster membership
    # -----------------------------------------------------------------------

    def get_cluster_members(
        self,
        uniref90_id: str,
        fetch_sequences: bool = False,
    ) -> pd.DataFrame:
        """
        List all member accessions of a UniRef90 cluster.

        Parameters
        ----------
        fetch_sequences : bool
            If True, retrieve the sequence for every member (may be slow —
            triggers UniProt API calls for unreviewed entries not yet cached).

        Returns
        -------
        DataFrame: member_accession, is_representative, sequence (optional),
                   interpro_ids, pfam_ids
        """
        with self._get_conn() as conn:
            rows = conn.execute(
                """SELECT cm.member_accession, cm.is_representative
                   FROM cluster_members cm
                   WHERE cm.uniref90_id = ?""",
                (uniref90_id,),
            ).fetchall()

        if not rows:
            return pd.DataFrame(columns=["member_accession", "is_representative"])

        df = pd.DataFrame([dict(r) for r in rows])

        if fetch_sequences:
            df["sequence"] = df["member_accession"].apply(self.fetch_member_sequence)

        # Attach domain summaries
        with self._get_conn() as conn:
            accs = df["member_accession"].tolist()
            ph   = ",".join("?" * len(accs))
            dom_rows = conn.execute(
                f"""SELECT uniprot_accession,
                           GROUP_CONCAT(DISTINCT interpro_id)  AS interpro_ids,
                           GROUP_CONCAT(DISTINCT CASE WHEN db_name='Pfam'
                                        THEN db_accession END) AS pfam_ids
                    FROM domains WHERE uniprot_accession IN ({ph})
                    GROUP BY uniprot_accession""",
                accs,
            ).fetchall()
        dom_df = pd.DataFrame([dict(r) for r in dom_rows]).rename(
            columns={"uniprot_accession": "member_accession"}
        )
        return df.merge(dom_df, on="member_accession", how="left").reset_index(drop=True)

    def get_cluster_for_accession(self, uniprot_acc: str) -> Dict:
        """
        Return the UniRef90 cluster that contains a given UniProtKB accession.

        Returns
        -------
        dict: uniref90_id, rep_accession, member_count, is_representative
        """
        with self._get_conn() as conn:
            id_row = conn.execute(
                "SELECT uniref90_id FROM idmapping WHERE uniprot_accession = ? LIMIT 1",
                (uniprot_acc,),
            ).fetchone()
            if not id_row:
                return {}

            uniref90_id = id_row["uniref90_id"]

            seq_row = conn.execute(
                "SELECT rep_accession, member_count FROM sequences WHERE uniref90_id = ?",
                (uniref90_id,),
            ).fetchone()

            mem_row = conn.execute(
                """SELECT is_representative FROM cluster_members
                   WHERE uniref90_id = ? AND member_accession = ?""",
                (uniref90_id, uniprot_acc),
            ).fetchone()

        return {
            "uniref90_id":       uniref90_id,
            "rep_accession":     seq_row["rep_accession"] if seq_row else None,
            "member_count":      seq_row["member_count"]  if seq_row else None,
            "is_representative": bool(mem_row["is_representative"]) if mem_row else False,
        }

    # -----------------------------------------------------------------------
    # On-demand sequence fetching
    # -----------------------------------------------------------------------

    def fetch_member_sequence(self, uniprot_accession: str) -> Optional[str]:
        """
        Return the amino acid sequence for a UniProtKB accession.

        Checks sequences_cache first; fetches from the UniProt REST API on
        a miss and caches the result.  Returns None gracefully on 404.

        Rate-limited to ≤3 requests/second (sleeps 0.34 s between calls).
        """
        with self._get_conn() as conn:
            row = conn.execute(
                "SELECT sequence FROM sequences_cache WHERE accession = ?",
                (uniprot_accession,),
            ).fetchone()
            if row:
                return row["sequence"]

        # Fetch from UniProt REST API
        url = _UNIPROT_FASTA_URL.format(accession=uniprot_accession)
        try:
            resp = requests.get(url, timeout=30)
            time.sleep(_RATE_LIMIT_SLEEP)
        except requests.exceptions.RequestException as exc:
            log.warning("Network error fetching %s: %s", uniprot_accession, exc)
            return None

        if resp.status_code == 404:
            return None
        resp.raise_for_status()

        fasta = resp.text
        seq = "".join(
            line for line in fasta.splitlines() if line and not line.startswith(">")
        )
        if not seq:
            return None

        today = date.today().isoformat()
        with self._get_conn() as conn:
            with conn:
                conn.execute(
                    """INSERT OR REPLACE INTO sequences_cache
                       (accession, sequence, reviewed, fetched_date)
                       VALUES (?, ?, 0, ?)""",
                    (uniprot_accession, seq, today),
                )
        return seq
