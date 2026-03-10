#!/usr/bin/env bash
# load_to_shm.sh — Copy/decompress search index files into fast local storage.
#
# Run once after initial setup, or if files are lost/corrupted.
# Files are stored on local disk (/opt/shared/jhilzinger/uniprot/) and
# persist across reboots — no need to re-run after reboot.
#
# Usage: bash load_to_shm.sh [--force]
#   --force : overwrite existing files even if they already exist
#
# Files placed in /opt/shared/jhilzinger/uniprot/:
#   uniref90_diamond.dmnd  (~86 GB)   — copied from NFS diamond_db/
#   uniref50.fasta         (~38 GB)   — decompressed from NFS data/uniref50.fasta.gz
# Total footprint: ~124 GB

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FAST_DIR="/opt/shared/jhilzinger/uniprot"
LOG_FILE="$SCRIPT_DIR/build.log"
DIAMOND_SRC="$SCRIPT_DIR/diamond_db/uniref90_diamond.dmnd"
UNIREF50_GZ="$SCRIPT_DIR/data/uniref50.fasta.gz"
DIAMOND_DST="$FAST_DIR/uniref90_diamond.dmnd"
UNIREF50_DST="$FAST_DIR/uniref50.fasta"
FORCE=0
REQUIRED_GB=130   # ~124 GB needed + 6 GB headroom

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
for arg in "$@"; do
    [[ "$arg" == "--force" ]] && FORCE=1
done

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
log "load_to_shm.sh: starting — target: $FAST_DIR"

if [[ ! -f "$DIAMOND_SRC" ]]; then
    log "ERROR: DIAMOND index not found: $DIAMOND_SRC"
    exit 1
fi

if [[ ! -f "$UNIREF50_GZ" ]]; then
    log "ERROR: UniRef50 FASTA not found: $UNIREF50_GZ"
    log "       Download with:"
    log "       wget -c -O $UNIREF50_GZ https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
    exit 1
fi

# Check available space on target filesystem
FREE_KB=$(df -k "$FAST_DIR" 2>/dev/null | awk 'NR==2{print $4}' || df -k /opt/shared | awk 'NR==2{print $4}')
FREE_GB=$(( FREE_KB / 1024 / 1024 ))
if [[ "$FREE_GB" -lt "$REQUIRED_GB" ]]; then
    log "ERROR: insufficient space on $(df -k /opt/shared | awk 'NR==2{print $1}'). Need ~${REQUIRED_GB} GB, have ${FREE_GB} GB free."
    exit 1
fi
log "Preflight OK — ${FREE_GB} GB free on target filesystem (need ~${REQUIRED_GB} GB)"

# Create target directory
mkdir -p "$FAST_DIR"

T0=$(date +%s)

# ---------------------------------------------------------------------------
# Copy DIAMOND index (~86 GB)
# ---------------------------------------------------------------------------
if [[ -f "$DIAMOND_DST" && -s "$DIAMOND_DST" && "$FORCE" -eq 0 ]]; then
    EXISTING_GB=$(( $(stat -c%s "$DIAMOND_DST") / 1024 / 1024 / 1024 ))
    log "DIAMOND index already present (${EXISTING_GB} GB) — skipping (use --force to overwrite)"
else
    log "Copying DIAMOND index (~86 GB) ..."
    T1=$(date +%s)
    cp "$DIAMOND_SRC" "$DIAMOND_DST"
    T2=$(date +%s)
    SIZE_GB=$(( $(stat -c%s "$DIAMOND_DST") / 1024 / 1024 / 1024 ))
    log "DIAMOND index copied: ${SIZE_GB} GB in $(( T2 - T1 )) s"
fi

# ---------------------------------------------------------------------------
# Decompress UniRef50 FASTA (~38 GB uncompressed)
# ---------------------------------------------------------------------------
if [[ -f "$UNIREF50_DST" && -s "$UNIREF50_DST" && "$FORCE" -eq 0 ]]; then
    EXISTING_GB=$(( $(stat -c%s "$UNIREF50_DST") / 1024 / 1024 / 1024 ))
    log "UniRef50 FASTA already present (${EXISTING_GB} GB) — skipping (use --force to overwrite)"
else
    log "Decompressing UniRef50 FASTA (~38 GB uncompressed) ..."
    T1=$(date +%s)
    gunzip -c "$UNIREF50_GZ" > "$UNIREF50_DST"
    T2=$(date +%s)
    SIZE_GB=$(( $(stat -c%s "$UNIREF50_DST") / 1024 / 1024 / 1024 ))
    log "UniRef50 FASTA decompressed: ${SIZE_GB} GB in $(( T2 - T1 )) s"
fi

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
ERRORS=0
for f in "$DIAMOND_DST" "$UNIREF50_DST"; do
    if [[ ! -f "$f" || ! -s "$f" ]]; then
        log "ERROR: missing or empty: $f"
        ERRORS=$(( ERRORS + 1 ))
    fi
done
if [[ "$ERRORS" -gt 0 ]]; then
    log "Validation FAILED — $ERRORS file(s) missing"
    exit 1
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
T_END=$(date +%s)
WALL=$(( T_END - T0 ))
log "load_to_shm.sh: complete in ${WALL}s ($(( WALL / 60 )) min)"
log "$(df -h /opt/shared | tail -1)"
ls -lh "$FAST_DIR" | tee -a "$LOG_FILE"
