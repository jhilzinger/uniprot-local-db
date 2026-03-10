# jackhmmer NFS Benchmark — thar server (2026-03-05)

## Setup
- Server: thar, NFS-mounted storage at `/auto/sahara`
- jackhmmer: HMMER 3.3.2 at `/usr/bin/jackhmmer`
- Database: Swiss-Prot extracted from `uniprot_sprot.dat.gz`
  - 574,627 proteins, **0.21 GB** (`data/sprot_test.fasta`)
- Query: P04637 (human TP53, 117 aa N-terminal fragment)

## Wall Times

| Test | Config | Wall time |
|------|--------|-----------|
| A | 1 iteration, 8 threads | 19.2 s |
| B | 3 iterations, 8 threads (production default) | 4.0 s |
| C | 3 iterations, 1 thread (I/O baseline) | 9.6 s |

Note: Tests B and C ran only 2 effective iterations (no new seeds found on iteration 3).

## Hits Found

| Test | Hits (E-value ≤ 0.001) |
|------|------------------------|
| A (1 iter) | 34 |
| B (3 iter, 8 cpu) | 54 |
| C (3 iter, 1 cpu) | 54 |

Top hits for P04637: P04637, Q9TTA1, P56423, P56424, P61260 — all TP53/p53 family. ✓

## NFS I/O Rate
- Test B (8 threads): 204 MB / 4.0 s ≈ **51 MB/sec**
- Test C (1 thread):  204 MB / 9.6 s ≈ **21 MB/sec**
- Healthy NFS: 100–200 MB/sec. **NFS throughput is degraded (~21–51 MB/sec).**

## UniRef90 Projection

| Item | Value |
|------|-------|
| Swiss-Prot size | 0.21 GB |
| UniRef90 FASTA size | 84 GB |
| Scale factor | ~400× |
| Test B wall time | 4.0 s |
| **Projected UniRef90 wall time** | **~27 min** |

**Verdict: IMPRACTICAL** for routine use. Actual time likely worse — UniRef90's larger sequence space will trigger more iteration rounds than Swiss-Prot did.

## Key Findings
- jackhmmer is I/O-bound on NFS; 8 threads gives only ~2.4× speedup over 1 thread
- DIAMOND and MMseqs2 have the same underlying NFS problem (mmap-based index access)
- `/dev/shm` is available: 498 GB tmpfs, 498 GB free — could host DIAMOND DB (92 GB) + UniRef90 FASTA (84 GB) in RAM to bypass NFS entirely
- Copying to `/dev/shm` is a ~30 min one-time cost per reboot; contents lost on reboot

## Resolution (2026-03-09)

`/dev/shm` was rejected because its contents are lost on reboot. Instead, files are
stored on `/opt/shared/jhilzinger/uniprot/` (local XFS disk, 8.1 TB, sanctioned by
QB3 sysadmin), which persists across reboots. `load_fast_storage.sh` does the one-time
copy (~16 min). The jackhmmer target was also switched from UniRef90 (84 GB) to
UniRef50 (23 GB uncompressed) — 3.6× smaller, negligible sensitivity loss for twilight
zone detection. Actual measured runtimes: DIAMOND ~9.4 min, jackhmmer ~11.6 min (1 iter,
P04637 query, validated 2026-03-09).
