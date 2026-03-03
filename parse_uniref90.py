#!/usr/bin/env python3
"""
Parse UniRef90 FASTA (gzipped) into the SQLite `sequences` table.

Header format:
  >UniRef90_XXXXX cluster_name n=<count> Tax=<taxon> TaxID=<id> RepID=<rep>

Usage:
  python3 parse_uniref90.py <uniref90.fasta.gz> <db_path>
"""

import gzip
import logging
import re
import sqlite3
import sys
from typing import List, Tuple

from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BATCH_SIZE = 10_000

# Matches: >UniRef90_ID  cluster name  n=123  Tax=Org name  TaxID=9606  RepID=REPID
HEADER_RE = re.compile(
    r"^>(\S+)\s+(.*?)\s+n=(\d+)\s+Tax=.+?\s+TaxID=(\d+)\s+RepID=(\S+)$"
)


def open_db(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA cache_size=-2000000")   # 2 GB
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    return conn


def flush(conn: sqlite3.Connection, batch: List[Tuple]) -> None:
    with conn:
        conn.executemany(
            """INSERT OR REPLACE INTO sequences
               (uniref90_id, rep_accession, cluster_name, member_count, taxonomy_id, sequence)
               VALUES (?, ?, ?, ?, ?, ?)""",
            batch,
        )


def parse_and_insert(fasta_gz: str, db_path: str) -> None:
    conn = open_db(db_path)

    batch: List[Tuple] = []
    total = 0
    skipped = 0

    current_header: str | None = None
    seq_parts: List[str] = []

    log.info("Streaming %s …", fasta_gz)

    with gzip.open(fasta_gz, "rt") as fh:
        for raw_line in tqdm(fh, desc="UniRef90 lines", unit="lines", mininterval=10.0):
            line = raw_line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                # Emit previous record
                if current_header is not None:
                    m = HEADER_RE.match(current_header)
                    if m:
                        uid, name, n, taxid, repid = m.groups()
                        batch.append((uid, repid, name, int(n), int(taxid), "".join(seq_parts)))
                        total += 1
                        if len(batch) >= BATCH_SIZE:
                            flush(conn, batch)
                            batch.clear()
                            if total % 500_000 == 0:
                                log.info("  … %d sequences inserted", total)
                    else:
                        skipped += 1
                current_header = line
                seq_parts = []
            else:
                seq_parts.append(line)

    # Final record
    if current_header is not None:
        m = HEADER_RE.match(current_header)
        if m:
            uid, name, n, taxid, repid = m.groups()
            batch.append((uid, repid, name, int(n), int(taxid), "".join(seq_parts)))
            total += 1

    if batch:
        flush(conn, batch)

    conn.close()
    log.info("Done. Inserted %d sequences (%d headers skipped).", total, skipped)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <uniref90.fasta.gz> <db_path>")
    parse_and_insert(sys.argv[1], sys.argv[2])
