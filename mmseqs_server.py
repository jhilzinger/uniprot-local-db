#!/usr/bin/env python3
"""
MMseqs2 HTTP search server.

Accepts POST /search with a FASTA body and optional query-string parameters,
runs mmseqs easy-search against the pre-built UniRef90 database, and returns
tab-separated results (m8 format).

Compatible with the uniprotdb.py _search_via_server() client.

Usage:
    python3 mmseqs_server.py            # foreground
    nohup python3 mmseqs_server.py &    # background

Query parameters (all optional, defaults match uniprotdb.py):
    sensitivity   float   (default 7.5)
    max-seqs      int     (default 100)
    min-seq-id    float   (default 0.3)
    coverage      float   (default 0.5)
    format-output str     (default: query,target,pident,… m8 columns)
"""

import logging
import os
import shutil
import subprocess
import tempfile
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import parse_qs, urlparse

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
MMSEQS_DB  = os.path.join(BASE_DIR, "mmseqs_db", "uniref90_mmseqs")
PORT       = 8080
THREADS    = 8          # mmseqs threads per search
MAX_CONCURRENT = 4      # semaphore: cap parallel searches at 4 × THREADS

_DEFAULT_FORMAT = (
    "query,target,pident,alnlen,mismatch,gapopen,"
    "qstart,qend,tstart,tend,evalue,bits"
)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

_sem = threading.Semaphore(MAX_CONCURRENT)

# ---------------------------------------------------------------------------
# Search worker
# ---------------------------------------------------------------------------

def run_search(fasta: str, sensitivity: str, max_seqs: str,
               min_seq_id: str, coverage: str, fmt: str) -> str:
    tmp_dir = tempfile.mkdtemp(prefix="mmseqs_srv_")
    try:
        query_file  = os.path.join(tmp_dir, "query.fasta")
        result_file = os.path.join(tmp_dir, "result.tsv")
        work_dir    = os.path.join(tmp_dir, "work")
        os.makedirs(work_dir)

        with open(query_file, "w") as fh:
            fh.write(fasta)

        cmd = [
            "mmseqs", "easy-search",
            query_file, MMSEQS_DB, result_file, work_dir,
            "--format-output", fmt,
            "-s",           sensitivity,
            "--max-seqs",   max_seqs,
            "--min-seq-id", min_seq_id,
            "-c",           coverage,
            "--threads",    str(THREADS),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr[-2000:])
        if proc.stderr:
            log.debug("mmseqs stderr: %s", proc.stderr[:500])

        if os.path.exists(result_file):
            with open(result_file) as fh:
                return fh.read()
        return ""
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# HTTP handler
# ---------------------------------------------------------------------------

class Handler(BaseHTTPRequestHandler):

    def do_GET(self):
        if urlparse(self.path).path == "/health":
            self._respond(200, b"OK\n", "text/plain")
        else:
            self.send_error(404)

    def do_POST(self):
        parsed = urlparse(self.path)
        if parsed.path != "/search":
            self.send_error(404)
            return

        params = parse_qs(parsed.query)

        length = int(self.headers.get("Content-Length", 0))
        fasta  = self.rfile.read(length).decode("utf-8", errors="replace")

        if not fasta.strip():
            self.send_error(400, "Empty FASTA body")
            return

        sensitivity = params.get("sensitivity", ["7.5"])[0]
        max_seqs    = params.get("max-seqs",    ["100"])[0]
        min_seq_id  = params.get("min-seq-id",  ["0.3"])[0]
        coverage    = params.get("coverage",    ["0.5"])[0]
        fmt         = params.get("format-output", [_DEFAULT_FORMAT])[0]

        log.info("search request: s=%s max-seqs=%s min-id=%s cov=%s",
                 sensitivity, max_seqs, min_seq_id, coverage)

        try:
            with _sem:
                result = run_search(fasta, sensitivity, max_seqs,
                                    min_seq_id, coverage, fmt)
        except Exception as exc:
            log.error("search failed: %s", exc)
            self.send_error(500, str(exc)[:200])
            return

        body = result.encode("utf-8")
        self._respond(200, body, "text/plain; charset=utf-8")

    def _respond(self, code: int, body: bytes, content_type: str):
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt, *args):  # suppress default access log spam
        pass


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if not os.path.exists(MMSEQS_DB + ".dbtype"):
        raise SystemExit(
            f"ERROR: MMseqs2 database not found at {MMSEQS_DB}\n"
            "       Run build_database.sh first."
        )

    log.info("Starting MMseqs2 search server on port %d", PORT)
    log.info("  Database : %s", MMSEQS_DB)
    log.info("  Threads per search : %d  (max concurrent: %d)", THREADS, MAX_CONCURRENT)
    log.info("  Endpoint : http://localhost:%d/search", PORT)
    log.info("  Health   : http://localhost:%d/health", PORT)

    server = HTTPServer(("localhost", PORT), Handler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        log.info("Shutting down.")
