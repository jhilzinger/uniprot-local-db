#!/usr/bin/env bash
# =============================================================================
# build_database.sh — end-to-end UniProt local reference database build
#
# Checkpointed: each step is skipped if its .step_complete flag exists.
# Resume a partial build by simply re-running this script.
#
# Usage:
#   bash build_database.sh
# =============================================================================

set -euo pipefail

# Ensure conda-managed binaries (mmseqs, aria2c, python) are on PATH
# regardless of whether the invoking shell activated the conda environment.
export PATH="/usr2/people/jhilzinger/miniforge3/bin:$PATH"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR="/auto/sahara/namib/home/jhilzinger/databases/uniprot_db"
DATA_DIR="$BASE_DIR/data"
MMSEQS_DIR="$BASE_DIR/mmseqs_db"
DB_PATH="$BASE_DIR/uniprot.sqlite"
LOG="$BASE_DIR/build.log"

CONDA="/usr2/people/jhilzinger/miniforge3/bin/conda"
PYTHON="/usr2/people/jhilzinger/miniforge3/bin/python3"
PIP="/usr2/people/jhilzinger/miniforge3/bin/pip"

START_TIME=$(date +%s)
mkdir -p "$DATA_DIR" "$MMSEQS_DIR"

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOG"
}

log_disk() {
    local label="$1"
    log "--- Disk usage ($label) ---"
    df -h "$BASE_DIR" | tee -a "$LOG"
    du -sh "$DATA_DIR" "$MMSEQS_DIR" "$DB_PATH" 2>/dev/null | tee -a "$LOG" || true
}

check_step()  { [ -f "$BASE_DIR/.$1_complete" ]; }
mark_step()   { touch "$BASE_DIR/.$1_complete"; log "✓ Step '$1' complete."; }

# ---------------------------------------------------------------------------
# Step 0 — Dependency checks
# ---------------------------------------------------------------------------
log "=== Checking dependencies ==="

# MMseqs2
if ! command -v mmseqs &>/dev/null; then
    log "MMseqs2 not found — installing via conda …"
    "$CONDA" install -y -c conda-forge -c bioconda mmseqs2 2>&1 | tee -a "$LOG"
fi
MMSEQS_VER=$(mmseqs version 2>&1 | head -1)
log "  mmseqs  : $MMSEQS_VER"

# aria2c
if ! command -v aria2c &>/dev/null; then
    log "aria2c not found — installing via conda …"
    "$CONDA" install -y -c conda-forge aria2 2>&1 | tee -a "$LOG"
fi
ARIA_VER=$(aria2c --version 2>&1 | head -1)
log "  aria2c  : $ARIA_VER"

# Python packages
log "Ensuring Python packages: biopython pandas requests tqdm …"
"$PIP" install --quiet biopython pandas requests tqdm 2>&1 | tee -a "$LOG"
PY_VER=$("$PYTHON" --version 2>&1)
log "  python  : $PY_VER"

log "All dependencies satisfied."
log_disk "before download"

# ---------------------------------------------------------------------------
# Step 1 — Download
# ---------------------------------------------------------------------------
if ! check_step download; then
    log "=== Step 1: Download source files ==="
    cd "$DATA_DIR"

    # UniProt no longer ships .sum files; MD5s are embedded in RELEASE.metalink XML.
    # Download metalink XML files with curl (NOT aria2c, which auto-processes them
    # as download manifests and would pull the entire release directory).
    declare -A METALINKS=(
        ["uniref90.RELEASE.metalink"]="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/RELEASE.metalink"
        ["sprot.RELEASE.metalink"]="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/RELEASE.metalink"
        ["idmapping.RELEASE.metalink"]="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/RELEASE.metalink"
    )
    for filename in "${!METALINKS[@]}"; do
        url="${METALINKS[$filename]}"
        if [ -f "$filename" ]; then
            log "  Already present: $filename — skipping."
        else
            log "  Downloading (curl): $filename"
            curl -sSL -o "$filename" "$url" 2>&1 | tee -a "$LOG"
        fi
    done

    # Large data files — use aria2c for resume support and parallel connections.
    declare -A FILES=(
        ["uniref90.fasta.gz"]="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
        ["uniprot_sprot.dat.gz"]="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"
        ["idmapping_selected.tab.gz"]="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"
        ["protein2ipr.dat.gz"]="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz"
    )
    for filename in "${!FILES[@]}"; do
        url="${FILES[$filename]}"
        if [ -f "$filename" ]; then
            log "  Already present: $filename — skipping."
        else
            log "  Downloading: $filename"
            aria2c -x4 -c \
                --log="$LOG" --log-level=notice \
                -o "$filename" "$url" 2>&1 | tee -a "$LOG"
        fi
    done

    mark_step download
    log_disk "after download"
fi

# ---------------------------------------------------------------------------
# Step 2 — Checksum verification (via RELEASE.metalink XML)
# ---------------------------------------------------------------------------
if ! check_step checksum; then
    log "=== Step 2: Verify checksums ==="

    "$PYTHON" - <<PYEOF
import xml.etree.ElementTree as ET, hashlib, sys, os

DATA_DIR = "$DATA_DIR"
LOG      = "$LOG"

# Map local filename → (metalink_file, name_in_metalink)
CHECKS = [
    ("uniref90.fasta.gz",        "uniref90.RELEASE.metalink",    "uniref90.fasta.gz"),
    ("uniprot_sprot.dat.gz",     "sprot.RELEASE.metalink",       "uniprot_sprot.dat.gz"),
    ("idmapping_selected.tab.gz","idmapping.RELEASE.metalink",   "idmapping_selected.tab.gz"),
]

NS = "http://www.metalinker.org/"

def extract_md5(metalink_path, target_name):
    tree = ET.parse(metalink_path)
    root = tree.getroot()
    for f in root.iter(f"{{{NS}}}file"):
        if f.get("name") == target_name:
            for h in f.iter(f"{{{NS}}}hash"):
                if h.get("type") == "md5":
                    return h.text.strip().lower()
    return None

def md5_file(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()

fail = False
for local_name, metalink_name, entry_name in CHECKS:
    local_path    = os.path.join(DATA_DIR, local_name)
    metalink_path = os.path.join(DATA_DIR, metalink_name)

    if not os.path.exists(local_path):
        print(f"  SKIP (not downloaded yet): {local_name}")
        continue
    if not os.path.exists(metalink_path):
        print(f"  SKIP (no metalink): {local_name}")
        continue

    expected = extract_md5(metalink_path, entry_name)
    if not expected:
        print(f"  WARNING: MD5 not found in metalink for {entry_name}")
        continue

    print(f"  Verifying {local_name} …", flush=True)
    actual = md5_file(local_path)
    if actual == expected:
        msg = f"  OK: {local_name}  ({actual})"
    else:
        msg = f"  MISMATCH: {local_name}\n    expected {expected}\n    got      {actual}"
        fail = True
    print(msg)
    with open(LOG, "a") as lf:
        lf.write(msg + "\n")

# protein2ipr has no metalink — log size only
p2ipr = os.path.join(DATA_DIR, "protein2ipr.dat.gz")
if os.path.exists(p2ipr):
    size = os.path.getsize(p2ipr)
    print(f"  protein2ipr.dat.gz : {size:,} bytes (no checksum available)")

if fail:
    print("ERROR: Checksum mismatch — delete corrupted file(s) and re-run.")
    sys.exit(1)
PYEOF

    mark_step checksum
fi

# ---------------------------------------------------------------------------
# Step 3 — SQLite initialisation
# ---------------------------------------------------------------------------
if ! check_step sqlite_init; then
    log "=== Step 3: Initialise SQLite database ==="

    "$PYTHON" - <<'PYEOF'
import sqlite3, sys

db_path = "/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/uniprot.sqlite"
conn = sqlite3.connect(db_path)
conn.execute("PRAGMA journal_mode=WAL")
conn.execute("PRAGMA cache_size=-2000000")  # 2 GB

conn.executescript("""
CREATE TABLE IF NOT EXISTS sequences (
    uniref90_id     TEXT PRIMARY KEY,
    rep_accession   TEXT NOT NULL,
    cluster_name    TEXT,
    member_count    INT,
    taxonomy_id     INT,
    sequence        TEXT
);

CREATE TABLE IF NOT EXISTS cluster_members (
    uniref90_id        TEXT NOT NULL,
    member_accession   TEXT NOT NULL,
    is_representative  BOOLEAN NOT NULL DEFAULT 0,
    PRIMARY KEY (uniref90_id, member_accession)
);

CREATE TABLE IF NOT EXISTS idmapping (
    uniprot_accession  TEXT NOT NULL,
    uniref90_id        TEXT,
    uniref50_id        TEXT,
    uniref100_id       TEXT
);

CREATE TABLE IF NOT EXISTS domains (
    uniprot_accession  TEXT NOT NULL,
    interpro_id        TEXT,
    interpro_name      TEXT,
    db_name            TEXT,
    db_accession       TEXT,
    start_pos          INT,
    end_pos            INT
);

CREATE TABLE IF NOT EXISTS sprot_features (
    uniprot_accession  TEXT NOT NULL,
    feature_type       TEXT NOT NULL,
    start_pos          INT,
    end_pos            INT,
    description        TEXT
);

CREATE TABLE IF NOT EXISTS sequences_cache (
    accession    TEXT PRIMARY KEY,
    sequence     TEXT NOT NULL,
    reviewed     BOOLEAN,
    fetched_date TEXT
);

-- Indexes
CREATE INDEX IF NOT EXISTS idx_sequences_rep         ON sequences(rep_accession);
CREATE INDEX IF NOT EXISTS idx_cluster_members_uniref    ON cluster_members(uniref90_id);
CREATE INDEX IF NOT EXISTS idx_cluster_members_accession ON cluster_members(member_accession);
CREATE INDEX IF NOT EXISTS idx_idmapping_accession    ON idmapping(uniprot_accession);
CREATE INDEX IF NOT EXISTS idx_idmapping_uniref90     ON idmapping(uniref90_id);
CREATE INDEX IF NOT EXISTS idx_domains_accession      ON domains(uniprot_accession);
CREATE INDEX IF NOT EXISTS idx_domains_interpro       ON domains(interpro_id);
CREATE INDEX IF NOT EXISTS idx_domains_db             ON domains(db_name, db_accession);
CREATE INDEX IF NOT EXISTS idx_sprot_features_accession ON sprot_features(uniprot_accession);
CREATE INDEX IF NOT EXISTS idx_sprot_features_type    ON sprot_features(feature_type);
""")

conn.commit()
conn.close()
print("SQLite database initialised successfully.")
PYEOF

    mark_step sqlite_init
fi

# ---------------------------------------------------------------------------
# Step 4 — Parse UniRef90 FASTA
# ---------------------------------------------------------------------------
if ! check_step parse_uniref90; then
    log "=== Step 4: Parse UniRef90 FASTA ==="
    "$PYTHON" "$BASE_DIR/parse_uniref90.py" \
        "$DATA_DIR/uniref90.fasta.gz"       \
        "$DB_PATH"                          \
        2>&1 | tee -a "$LOG"
    mark_step parse_uniref90
    log_disk "after parse_uniref90"
fi

# ---------------------------------------------------------------------------
# Step 5 — Parse ID mapping
# ---------------------------------------------------------------------------
if ! check_step parse_idmapping; then
    log "=== Step 5: Parse ID mapping ==="
    "$PYTHON" "$BASE_DIR/parse_idmapping.py"       \
        "$DATA_DIR/idmapping_selected.tab.gz"      \
        "$DB_PATH"                                 \
        2>&1 | tee -a "$LOG"
    mark_step parse_idmapping
    log_disk "after parse_idmapping"
fi

# ---------------------------------------------------------------------------
# Step 6 — Parse InterPro domain annotations
# ---------------------------------------------------------------------------
if ! check_step parse_protein2ipr; then
    log "=== Step 6: Parse InterPro domain annotations ==="
    "$PYTHON" "$BASE_DIR/parse_protein2ipr.py" \
        "$DATA_DIR/protein2ipr.dat.gz"         \
        "$DB_PATH"                             \
        2>&1 | tee -a "$LOG"
    mark_step parse_protein2ipr
    log_disk "after parse_protein2ipr"
fi

# ---------------------------------------------------------------------------
# Step 7 — Parse Swiss-Prot features
# ---------------------------------------------------------------------------
if ! check_step parse_sprot; then
    log "=== Step 7: Parse Swiss-Prot features ==="
    "$PYTHON" "$BASE_DIR/parse_sprot.py"      \
        "$DATA_DIR/uniprot_sprot.dat.gz"      \
        "$DB_PATH"                            \
        2>&1 | tee -a "$LOG"
    mark_step parse_sprot
    log_disk "after parse_sprot"
fi

# ---------------------------------------------------------------------------
# Step 8 — MMseqs2 database + index build
# ---------------------------------------------------------------------------
if ! check_step mmseqs_build; then
    log "=== Step 8: Build MMseqs2 database and index ==="
    log "  This step takes several hours; safe to interrupt and resume."

    TMP_MMSEQS="$MMSEQS_DIR/tmp"
    mkdir -p "$TMP_MMSEQS"

    log "  Creating MMseqs2 sequence database …"
    mmseqs createdb \
        "$DATA_DIR/uniref90.fasta.gz" \
        "$MMSEQS_DIR/uniref90_mmseqs" \
        2>&1 | tee -a "$LOG"

    log "  Creating MMseqs2 prefilter index (--split-memory-limit 500G, --threads 32) …"
    mmseqs createindex \
        "$MMSEQS_DIR/uniref90_mmseqs" \
        "$TMP_MMSEQS"                 \
        --threads 32                  \
        --split-memory-limit 500G     \
        --remove-tmp-files 1          \
        2>&1 | tee -a "$LOG"

    mark_step mmseqs_build
    log_disk "after mmseqs_build"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))
HOURS=$(( ELAPSED / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

log ""
log "================================================================"
log "  Build complete!"
log "  Elapsed: ${HOURS}h ${MINS}m ${SECS}s"
log ""
log "  Database : $DB_PATH"
log "  MMseqs2  : $MMSEQS_DIR/uniref90_mmseqs"
log ""
log "  Start the search server with:"
log "    nohup bash $BASE_DIR/mmseqs_server.sh &"
log "================================================================"
log_disk "final"
