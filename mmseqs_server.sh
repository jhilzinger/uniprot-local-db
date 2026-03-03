#!/usr/bin/env bash
# Start the MMseqs2 HTTP search server.
#
# Usage:
#   bash mmseqs_server.sh                      (foreground)
#   nohup bash mmseqs_server.sh &              (background, detached)
#
# The server exposes a REST endpoint at http://localhost:8080/search
# Query with: curl -s -X POST http://localhost:8080/search -d @query.fasta
# Health:     curl -s http://localhost:8080/health

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="/usr2/people/jhilzinger/miniforge3/bin:$PATH"

exec python3 "$BASE_DIR/mmseqs_server.py" "$@"
