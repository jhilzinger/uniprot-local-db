# UniProt Local Database — LLM Context

## Paths
- SQLite DB: `/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/uniprot.sqlite`
- Scripts: `/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/`
- Python API: `uniprotdb.py`
- Fast storage (local disk): `/opt/shared/jhilzinger/uniprot/` — populated once by `load_to_shm.sh`; persists across reboots

## Data (UniProt/InterPro release 2026-02-25)
- 188,848,220 UniRef90 cluster representatives with sequences
- 203,130,941 UniProtKB → UniRef90/50/100 ID mappings
- 3,532,456 Swiss-Prot functional feature annotations (reviewed entries only)
- 1,174,000,000 InterPro domain annotations (completed 2026-03-09)

## SQLite Schema

```sql
sequences        (uniref90_id PK, rep_accession, cluster_name, member_count, taxonomy_id, sequence)
cluster_members  (uniref90_id, member_accession, is_representative)
idmapping        (uniprot_accession, uniref90_id, uniref50_id, uniref100_id)
domains          (uniprot_accession, interpro_id, interpro_name, db_name, db_accession, start_pos, end_pos)
sprot_features   (uniprot_accession, feature_type, start_pos, end_pos, description)
sequences_cache  (accession PK, sequence, reviewed, fetched_date)
```

Key indexes: sequences.rep_accession, idmapping.uniprot_accession, idmapping.uniref90_id,
domains.uniprot_accession, domains.interpro_id, domains.db_name+db_accession

## Python API

```python
import sys
sys.path.insert(0, "/auto/sahara/namib/home/jhilzinger/databases/uniprot_db")
from uniprotdb import UniProtDB

db = UniProtDB(
    db_path="/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/uniprot.sqlite",
    # fast_storage_path defaults to /opt/shared/jhilzinger/uniprot — omit if using default
)

# Fast search — DIAMOND against UniRef90 in /dev/shm (~2–8 min)
hits = db.search_by_sequence(fasta_str, mode="fast", sensitivity="sensitive")
hits = db.search_by_sequence(fasta_str, mode="fast", sensitivity="fast", max_hits=10)
# → DataFrame: uniref90_id, rep_accession, evalue, identity, coverage,
#              interpro_ids, pfam_ids, member_count, taxonomy_id

# Sensitive search — jackhmmer against UniRef50 in /dev/shm (~5–15 min)
hits = db.search_by_sequence(fasta_str, mode="sensitive")
hits = db.search_by_sequence(fasta_str, mode="sensitive", jackhmmer_iterations=1, max_hits=10)
# → same DataFrame columns; identity/coverage are NaN
# UniRef50 hits are resolved to UniRef90 clusters via idmapping.uniref50_id

# Annotation queries (SQLite only, no /dev/shm required)
info  = db.get_by_accession("P04637")             # → dict: sequence, domains, sprot_features, ...
prots = db.search_by_domain(pfam_id="PF00069")    # → DataFrame
arch  = db.get_domain_architecture("P04637")      # → list of domain dicts
mems  = db.get_cluster_members("UniRef90_P04637") # → DataFrame
clust = db.get_cluster_for_accession("P04637")    # → dict
seq   = db.fetch_member_sequence("P04637")        # → str (cache + UniProt REST fallback)
status = db.search_status()                        # → dict; check shm_ready before searching
```

## Populate fast storage (run once)

```bash
bash /auto/sahara/namib/home/jhilzinger/databases/uniprot_db/load_to_shm.sh
# ~45 min wall time; ~124 GB copied/decompressed into /opt/shared/jhilzinger/uniprot/
# Files persist across reboots — only needs to be re-run if files are lost/corrupted
# Check readiness: db.search_status()["shm_ready"]
```

## Search architecture

| Mode | Engine | Target | Fast storage size | Notes |
|------|--------|--------|-------------------|-------|
| `mode="fast"` | DIAMOND | UniRef90 | ~86 GB (.dmnd) | identity/coverage populated |
| `mode="sensitive"` | jackhmmer | UniRef50 | ~38 GB (.fasta) | identity/coverage NaN |

**Why UniRef50 for jackhmmer:** jackhmmer is for twilight zone (<30% identity) where
50% vs 90% clustering threshold is irrelevant. UniRef50 is ~38 GB vs 84 GB (UniRef90),
2.2x smaller. UniRef50 hits map back to UniRef90 clusters via idmapping.uniref50_id.

**Fast storage:** /opt/shared/jhilzinger/uniprot/ (local XFS disk, 8.1 TB total,
2.9 TB free). ~124 GB footprint. Persists across reboots.

## NFS files (permanent, source of truth)

| File | Size | Notes |
|------|------|-------|
| `uniprot.sqlite` | ~331 GB | Annotation database |
| `diamond_db/uniref90_diamond.dmnd` | ~86 GB | DIAMOND index (NFS copy; runtime copy in /opt/shared) |
| `data/uniref90.fasta.gz` | ~45 GB | Compressed UniRef90 |
| `data/uniref50.fasta.gz` | ~12 GB | Compressed UniRef50 (source for /dev/shm decompression) |

## Known Issues / Status
- ✅ domains table: 1,174,000,000 rows, fully populated (completed 2026-03-09)
- ✅ sequences table: 188,848,220 rows
- ✅ sprot_features table: 3,532,456 rows
- ✅ idmapping table: 203,130,941 rows
- ✅ Fast storage confirmed: /opt/shared/jhilzinger/uniprot/ sanctioned by QB3 sysadmin 2026-03-09
- ✅ All annotation query methods functional (no fast storage required)
- Fast storage (/opt/shared/jhilzinger/uniprot/) persists across reboots — run load_to_shm.sh once
- sequences_cache: empty at build time, populated on demand via fetch_member_sequence()
- identity/coverage are NaN for jackhmmer results (not reported by --tblout)
- jackhmmer UniRef50 hits map back to UniRef90 via idmapping.uniref50_id (one extra join)

## Rebuild DIAMOND index

```bash
diamond makedb \
  --in data/uniref90.fasta.gz \
  --db diamond_db/uniref90_diamond \
  --threads 32
```

## Rebuild a single pipeline step

```bash
rm /auto/sahara/namib/home/jhilzinger/databases/uniprot_db/.parse_protein2ipr_complete
bash build_database.sh   # skips all other completed steps
```
