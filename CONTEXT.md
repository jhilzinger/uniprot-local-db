# UniProt Local Database — LLM Context

## Paths
- SQLite DB: `/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/uniprot.sqlite`
- DIAMOND db: `/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/diamond_db/uniref90_diamond.dmnd`
- Scripts: `/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/`
- Python API: `uniprotdb.py`

## Data (UniProt/InterPro release 2026-02-25)
- 188,848,220 UniRef90 cluster representatives with sequences
- 203,130,941 UniProtKB → UniRef90/50/100 ID mappings
- 3,532,456 Swiss-Prot functional feature annotations (reviewed proteins only)
- InterPro domain annotations: **1,174,000,000 rows** (re-parse completed 2026-03-09, bug fix applied 2026-03-03)

## SQLite Schema

```sql
sequences        (uniref90_id PK, rep_accession, cluster_name, member_count, taxonomy_id, sequence)
cluster_members  (uniref90_id, member_accession, is_representative)
idmapping        (uniprot_accession, uniref90_id, uniref50_id, uniref100_id)
domains          (uniprot_accession, interpro_id, interpro_name, db_name, db_accession, start_pos, end_pos)
sprot_features   (uniprot_accession, feature_type, start_pos, end_pos, description)
sequences_cache  (accession PK, sequence, reviewed, fetched_date)  -- populated on demand
```

Key indexes: sequences.rep_accession, idmapping.uniprot_accession, idmapping.uniref90_id,
domains.uniprot_accession, domains.interpro_id, domains.db_name+db_accession

## Python API

```python
from uniprotdb import UniProtDB
db = UniProtDB(
    db_path="/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/uniprot.sqlite",
    diamond_db_path="/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/diamond_db/uniref90_diamond.dmnd",
)

# Fast search (DIAMOND, 2–8 min depending on sensitivity)
hits   = db.search_by_sequence(fasta_str)                                           # mode="fast", sensitivity="sensitive" default
hits   = db.search_by_sequence(fasta_str, mode="fast", sensitivity="sensitive")     # → DataFrame: uniref90_id, rep_accession, evalue, identity, coverage, interpro_ids, pfam_ids, member_count, taxonomy_id
hits   = db.search_by_sequence(fasta_str, mode="fast", sensitivity="fast", max_hits=10)  # faster for testing

# Deep homology search (jackhmmer, iterative profile-profile, ~3-10 min)
hits   = db.search_by_sequence(fasta_str, mode="sensitive")                         # → same DataFrame columns; identity/coverage are NaN
hits   = db.search_by_sequence(fasta_str, mode="sensitive", jackhmmer_iterations=1, max_hits=10)  # faster for testing

info   = db.get_by_accession("P04637")             # → dict: sequence, domains, sprot_features, cluster_id, is_reviewed, cluster_member_count
prots  = db.search_by_domain(pfam_id="PF00069")    # → DataFrame: uniprot_accession, uniref90_id, domain_info, start_pos, end_pos, member_count
arch   = db.get_domain_architecture("P04637")      # → list of dicts: db_name, db_accession, interpro_id, name, start, end
mems   = db.get_cluster_members("UniRef90_P04637") # → DataFrame: member_accession, is_representative, interpro_ids, pfam_ids
clust  = db.get_cluster_for_accession("P04637")    # → dict: uniref90_id, rep_accession, member_count, is_representative
seq    = db.fetch_member_sequence("P04637")        # → str (checks cache, falls back to UniProt REST API)
status = db.diamond_status()                       # → dict: available, path, db_path, db_exists, db_size_gb, version, jackhmmer
```

UniRef90 uncompressed FASTA required for jackhmmer: `data/uniref90.fasta` (84 GB).

## DIAMOND Search
DIAMOND is used for fast homology search (replaces MMseqs2, which was NFS-incompatible).
No server process required — DIAMOND runs as a subprocess, NFS-safe.

Sensitivity modes (pass as `sensitivity=` parameter):

| Mode | Typical runtime | Use case |
|------|----------------|----------|
| `fast` | ~1–2 min | Quick screening |
| `sensitive` | ~3–5 min | Default, recommended |
| `more-sensitive` | ~5–8 min | Moderate divergence |
| `very-sensitive` | ~8–15 min | Distant homologs |
| `ultra-sensitive` | ~15–30 min | Approaches jackhmmer sensitivity |

Rebuild DIAMOND db if needed:
```bash
diamond makedb --in data/uniref90.fasta --db diamond_db/uniref90_diamond --threads 32
```

## Known Issues / Status
- `domains` table fully populated: 1,174,000,000 rows (completed 2026-03-09)
- `sequences_cache` empty at build time, populated on first `fetch_member_sequence()` / `get_by_accession()` call
- DIAMOND runs as subprocess (no persistent server); NFS-safe
- jackhmmer searches run synchronously; expect 3–10 min per query against UniRef90
- identity/coverage columns are NaN for jackhmmer results (not reported by --tblout)

## Rebuild a single step
```bash
rm /path/to/uniprot_db/.parse_protein2ipr_complete
bash build_database.sh   # skips all other completed steps
```
