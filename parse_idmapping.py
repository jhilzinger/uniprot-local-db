#!/usr/bin/env python3
"""
Parse idmapping_selected.tab.gz into the SQLite `idmapping` and
`cluster_members` tables.

Column layout (0-based):
  0  UniProtKB-AC      1  UniProtKB-ID     2  GeneID
  3  RefSeq            4  GI               5  PDB
  6  GO                7  UniRef100        8  UniRef90
  9  UniRef50          10 UniParc          11 PIR
  12 NCBI-taxon        ...

Usage:
  python3 parse_idmapping.py <idmapping_selected.tab.gz> <db_path>
"""

import gzip
import logging
import sqlite3
import sys
from typing import Dict, List, Tuple

from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BATCH_SIZE = 10_000

COL_AC       = 0
COL_REF100   = 7
COL_REF90    = 8
COL_REF50    = 9


def open_db(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA cache_size=-2000000")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    return conn


def load_cluster_rep_map(conn: sqlite3.Connection) -> Dict[str, str]:
    """
    Returns {uniref90_id: representative_accession}.

    The representative accession is derived from the cluster ID itself
    (UniRef90_P04637 → P04637), which is the canonical UniProtKB accession
    of the cluster representative.  This is more reliable than the RepID
    field in the FASTA header which may carry an entry name (ALBU_HUMAN).
    """
    log.info("Loading cluster IDs from sequences table …")
    result = {}
    for (uid,) in conn.execute("SELECT uniref90_id FROM sequences"):
        # "UniRef90_P04637" → "P04637"
        acc = uid.split("_", 1)[1]
        result[uid] = acc
    log.info("  Loaded %d cluster entries.", len(result))
    return result


def flush(
    conn: sqlite3.Connection,
    idmap_batch: List[Tuple],
    member_batch: List[Tuple],
) -> None:
    with conn:
        conn.executemany(
            """INSERT OR REPLACE INTO idmapping
               (uniprot_accession, uniref100_id, uniref90_id, uniref50_id)
               VALUES (?, ?, ?, ?)""",
            idmap_batch,
        )
        conn.executemany(
            """INSERT OR IGNORE INTO cluster_members
               (uniref90_id, member_accession, is_representative)
               VALUES (?, ?, ?)""",
            member_batch,
        )


def parse_and_insert(tab_gz: str, db_path: str) -> None:
    conn = open_db(db_path)
    cluster_rep = load_cluster_rep_map(conn)

    idmap_batch: List[Tuple] = []
    member_batch: List[Tuple] = []
    total = 0
    no_ref90 = 0

    log.info("Streaming %s …", tab_gz)

    with gzip.open(tab_gz, "rt") as fh:
        for raw_line in tqdm(fh, desc="idmapping lines", unit="lines", mininterval=10.0):
            cols = raw_line.rstrip("\n").split("\t")
            if len(cols) <= COL_REF50:
                continue

            ac       = cols[COL_AC]
            ref100   = cols[COL_REF100] or None
            ref90    = cols[COL_REF90]  or None
            ref50    = cols[COL_REF50]  or None

            idmap_batch.append((ac, ref100, ref90, ref50))

            if ref90:
                is_rep = 1 if cluster_rep.get(ref90) == ac else 0
                member_batch.append((ref90, ac, is_rep))
            else:
                no_ref90 += 1

            total += 1
            if total % BATCH_SIZE == 0:
                flush(conn, idmap_batch, member_batch)
                idmap_batch.clear()
                member_batch.clear()
                if total % 5_000_000 == 0:
                    log.info("  … %d rows processed", total)

    if idmap_batch or member_batch:
        flush(conn, idmap_batch, member_batch)

    conn.close()
    log.info(
        "Done. Processed %d rows (%d lacked UniRef90 mapping).",
        total, no_ref90,
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <idmapping_selected.tab.gz> <db_path>")
    parse_and_insert(sys.argv[1], sys.argv[2])
