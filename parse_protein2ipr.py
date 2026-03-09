#!/usr/bin/env python3
"""
Parse protein2ipr.dat.gz into the SQLite `domains` table.

Actual column layout (tab-separated, 0-based):
  0  UniProt accession
  1  InterPro ID       (e.g. IPR000001)
  2  InterPro name
  3  DB accession      (e.g. PF00155, TIGR01821, G3DSA:3.40.640.10)
  4  Start position    (1-based)
  5  End position      (1-based inclusive)

db_name is derived from the DB accession prefix.

Only rows whose UniProt accession is already present in the idmapping
table are inserted (limits table size to sequences we actually have
metadata for).

Usage:
  python3 parse_protein2ipr.py <protein2ipr.dat.gz> <db_path>
"""

import gzip
import logging
import sqlite3
import sys
from typing import List, Set, Tuple

from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BATCH_SIZE = 10_000

# Derive a human-readable db_name from the DB accession prefix.
# This is needed because uniprotdb.py queries: db_name = 'Pfam' AND db_accession = ?
_DB_PREFIXES = [
    ("PF",       "Pfam"),
    ("PTHR",     "PANTHER"),
    ("TIGR",     "TIGRFAM"),
    ("G3DSA:",   "CATH-Gene3D"),
    ("SSF",      "SUPERFAMILY"),
    ("PIRSF",    "PIRSF"),
    ("SM",       "SMART"),
    ("PS",       "PROSITE"),
    ("MF_",      "HAMAP"),
    ("SFLD",     "SFLD"),
    ("cd",       "CDD"),
    ("PR",       "PRINTS"),
]

def _db_name(db_acc: str) -> str:
    for prefix, name in _DB_PREFIXES:
        if db_acc.startswith(prefix):
            return name
    return db_acc.split(":")[0] if ":" in db_acc else db_acc


def open_db(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA cache_size=-2000000")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    return conn


def load_known_accessions(conn: sqlite3.Connection) -> Set[str]:
    log.info("Loading known UniProt accessions from idmapping table …")
    known: Set[str] = set()
    for (ac,) in conn.execute("SELECT uniprot_accession FROM idmapping"):
        known.add(ac)
    log.info("  %d accessions loaded.", len(known))
    return known


def flush(conn: sqlite3.Connection, batch: List[Tuple]) -> None:
    with conn:
        conn.executemany(
            """INSERT INTO domains
               (uniprot_accession, interpro_id, interpro_name,
                db_name, db_accession, start_pos, end_pos)
               VALUES (?, ?, ?, ?, ?, ?, ?)""",
            batch,
        )


def parse_and_insert(dat_gz: str, db_path: str) -> None:
    conn = open_db(db_path)
    known = load_known_accessions(conn)

    batch: List[Tuple] = []
    total = 0
    filtered = 0
    malformed = 0

    log.info("Streaming %s …", dat_gz)

    with gzip.open(dat_gz, "rt") as fh:
        for raw_line in tqdm(fh, desc="protein2ipr lines", unit="lines", mininterval=10.0):
            cols = raw_line.rstrip("\n").split("\t")
            if len(cols) < 6:
                malformed += 1
                continue

            ac = cols[0]
            if ac not in known:
                filtered += 1
                continue

            db_acc = cols[3] or None
            try:
                start = int(cols[4])
                end   = int(cols[5])
            except ValueError:
                malformed += 1
                continue

            batch.append((
                ac,
                cols[1] or None,          # interpro_id
                cols[2] or None,          # interpro_name
                _db_name(cols[3]),        # db_name  (derived from accession)
                db_acc,                   # db_accession
                start,
                end,
            ))
            total += 1

            if len(batch) >= BATCH_SIZE:
                flush(conn, batch)
                batch.clear()
                if total % 2_000_000 == 0:
                    log.info("  … %d domain rows inserted", total)

    if batch:
        flush(conn, batch)

    conn.close()
    log.info(
        "Done. Inserted %d domain rows (%d filtered out, %d malformed).",
        total, filtered, malformed,
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <protein2ipr.dat.gz> <db_path>")
    parse_and_insert(sys.argv[1], sys.argv[2])
