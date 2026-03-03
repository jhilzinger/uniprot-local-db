#!/usr/bin/env python3
"""
Parse uniprot_sprot.dat.gz (Swiss-Prot flat file) via BioPython into:
  - sprot_features  (structural / functional annotations)
  - sequences_cache (reviewed sequences, reviewed=1)

Feature types captured:
  SIGNAL, TRANSMEM, COILED, COMPBIAS, REGION, REPEAT, MOTIF,
  DOMAIN, PROPEP, TRANSIT, DISULFID, MOD_RES, BINDING, ACT_SITE,
  TOPO_DOM, ZN_FING, DNA_BIND

Usage:
  python3 parse_sprot.py <uniprot_sprot.dat.gz> <db_path>
"""

import gzip
import logging
import sqlite3
import sys
from datetime import date
from typing import List, Tuple

from Bio import SwissProt
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BATCH_SIZE = 10_000

FEATURE_TYPES = frozenset({
    "SIGNAL", "TRANSMEM", "COILED", "COMPBIAS", "REGION", "REPEAT",
    "MOTIF", "DOMAIN", "PROPEP", "TRANSIT", "DISULFID", "MOD_RES",
    "BINDING", "ACT_SITE", "TOPO_DOM", "ZN_FING", "DNA_BIND",
})

TODAY = date.today().isoformat()


def open_db(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA cache_size=-2000000")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    return conn


def flush(
    conn: sqlite3.Connection,
    feat_batch: List[Tuple],
    cache_batch: List[Tuple],
) -> None:
    with conn:
        conn.executemany(
            """INSERT OR IGNORE INTO sprot_features
               (uniprot_accession, feature_type, start_pos, end_pos, description)
               VALUES (?, ?, ?, ?, ?)""",
            feat_batch,
        )
        conn.executemany(
            """INSERT OR REPLACE INTO sequences_cache
               (accession, sequence, reviewed, fetched_date)
               VALUES (?, ?, 1, ?)""",
            cache_batch,
        )


def _get_description(qualifiers: dict) -> str:
    """Extract a human-readable description from a feature's qualifier dict."""
    for key in ("note", "description", "product"):
        val = qualifiers.get(key)
        if val:
            return val if isinstance(val, str) else str(val)
    return ""


def parse_and_insert(dat_gz: str, db_path: str) -> None:
    conn = open_db(db_path)

    feat_batch: List[Tuple] = []
    cache_batch: List[Tuple] = []
    rec_count = 0
    feat_count = 0

    log.info("Streaming %s …", dat_gz)

    with gzip.open(dat_gz, "rt") as fh:
        for record in tqdm(
            SwissProt.parse(fh),
            desc="Swiss-Prot records",
            unit="records",
            mininterval=10.0,
        ):
            ac = record.accessions[0]
            seq = record.sequence

            cache_batch.append((ac, seq, TODAY))

            for feat in record.features:
                ftype = feat.type
                if ftype not in FEATURE_TYPES:
                    continue

                # BioPython FeatureLocation is 0-based half-open; convert to
                # 1-based inclusive (UniProt convention)
                try:
                    start = int(feat.location.start) + 1
                    end   = int(feat.location.end)
                except (TypeError, AttributeError):
                    start = None
                    end   = None

                desc = _get_description(feat.qualifiers) if feat.qualifiers else ""
                feat_batch.append((ac, ftype, start, end, desc or None))
                feat_count += 1

            rec_count += 1
            if rec_count % BATCH_SIZE == 0:
                flush(conn, feat_batch, cache_batch)
                feat_batch.clear()
                cache_batch.clear()
                log.info("  … %d records, %d features", rec_count, feat_count)

    if feat_batch or cache_batch:
        flush(conn, feat_batch, cache_batch)

    conn.close()
    log.info(
        "Done. Parsed %d Swiss-Prot records, inserted %d feature rows.",
        rec_count, feat_count,
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <uniprot_sprot.dat.gz> <db_path>")
    parse_and_insert(sys.argv[1], sys.argv[2])
