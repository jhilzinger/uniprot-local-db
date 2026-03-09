# UniProt Local Database

A self-contained, locally-hosted reference database built from public UniProt/UniRef90
and InterPro data. Provides fast sequence similarity search and protein annotation lookup
without external network dependencies.

**Built:** 2026-03-03 | **Updated:** 2026-03-03 (DIAMOND migration) | **Runtime:** 137 h 12 m | **Total size:** ~323 GB on disk

---

## Contents

| Component | Description | Size on disk |
|-----------|-------------|-------------|
| `uniprot.sqlite` | Annotation database (SQLite) | 147 GB |
| `diamond_db/` | DIAMOND sequence index for fast search | 92 GB |
| `data/uniref90.fasta` | Uncompressed UniRef90 FASTA (required for jackhmmer) | 84 GB |

---

## Data Sources

| Source | File | Description |
|--------|------|-------------|
| **UniRef90** | `uniref90.fasta.gz` | Non-redundant protein clusters at 90% identity — 188.8 M sequences |
| **UniProt ID mapping** | `idmapping_selected.tab.gz` | Maps UniProtKB accessions to UniRef cluster IDs |
| **InterPro** | `protein2ipr.dat.gz` | Domain annotations (Pfam, PANTHER, CATH, etc.) — 1,174,000,000 rows |
| **Swiss-Prot** | `uniprot_sprot.dat.gz` | Curated feature annotations for reviewed entries |

All files are from the UniProt/InterPro release current as of 2026-02-25.

---

## SQLite Schema

### `sequences` — 188,848,220 rows
One row per UniRef90 cluster representative.

| Column | Type | Description |
|--------|------|-------------|
| `uniref90_id` | TEXT PK | UniRef90 cluster ID (e.g. `UniRef90_P04637`) |
| `rep_accession` | TEXT | UniProtKB accession of the cluster representative |
| `cluster_name` | TEXT | Cluster name from UniRef90 header |
| `member_count` | INT | Number of sequences in the cluster |
| `taxonomy_id` | INT | NCBI taxonomy ID of the representative |
| `sequence` | TEXT | Amino acid sequence of the representative |

### `cluster_members`
Maps individual UniProtKB accessions to their UniRef90 cluster.

| Column | Type | Description |
|--------|------|-------------|
| `uniref90_id` | TEXT | Foreign key → `sequences.uniref90_id` |
| `member_accession` | TEXT | UniProtKB accession |
| `is_representative` | BOOL | Whether this member is the cluster representative |

### `idmapping` — 203,130,941 rows
Cross-reference mapping between UniProtKB accessions and UniRef cluster IDs.

| Column | Type | Description |
|--------|------|-------------|
| `uniprot_accession` | TEXT | UniProtKB accession |
| `uniref90_id` | TEXT | Corresponding UniRef90 cluster |
| `uniref50_id` | TEXT | Corresponding UniRef50 cluster |
| `uniref100_id` | TEXT | Corresponding UniRef100 cluster |

### `domains` — 1,174,000,000 rows
Holds InterPro domain annotations per protein.

| Column | Type | Description |
|--------|------|-------------|
| `uniprot_accession` | TEXT | UniProtKB accession |
| `interpro_id` | TEXT | InterPro entry (e.g. `IPR000308`) |
| `interpro_name` | TEXT | Human-readable domain name |
| `db_name` | TEXT | Source database (Pfam, PANTHER, CATH-Gene3D, etc.) |
| `db_accession` | TEXT | Accession in the source database (e.g. `PF00076`) |
| `start_pos` | INT | Domain start position (1-based) |
| `end_pos` | INT | Domain end position |

### `sprot_features` — 3,532,456 rows
Functional site annotations from Swiss-Prot (reviewed entries only).

| Column | Type | Description |
|--------|------|-------------|
| `uniprot_accession` | TEXT | UniProtKB accession |
| `feature_type` | TEXT | Feature type (e.g. `Active site`, `Binding site`, `Disulfide bond`) |
| `start_pos` | INT | Feature start position |
| `end_pos` | INT | Feature end position |
| `description` | TEXT | Free-text annotation |

### `sequences_cache`
On-demand cache for sequences fetched from the UniProt REST API (populated at query time,
not during the build). Avoids re-fetching sequences already retrieved.

| Column | Type | Description |
|--------|------|-------------|
| `accession` | TEXT PK | UniProtKB accession |
| `sequence` | TEXT | Amino acid sequence |
| `reviewed` | BOOL | Whether the entry is Swiss-Prot reviewed |
| `fetched_date` | TEXT | ISO date when the sequence was fetched |

---

## Sequence Search

Two search engines are available, selectable via the `mode` parameter.

### Fast search — DIAMOND

Uses a precomputed DIAMOND protein database (`diamond_db/uniref90_diamond.dmnd`)
for searches that typically complete in 2–8 minutes depending on sensitivity mode.
DIAMOND performs sequential disk reads and is fully compatible with NFS storage.
No persistent server process required.

Sensitivity modes (pass as `sensitivity=` parameter):

| Mode | Typical runtime | Use case |
|------|----------------|----------|
| `fast` | ~1–2 min | Quick screening |
| `sensitive` | ~3–5 min | Default, recommended |
| `more-sensitive` | ~5–8 min | Moderate divergence |
| `very-sensitive` | ~8–15 min | Distant homologs |
| `ultra-sensitive` | ~15–30 min | Approaches jackhmmer sensitivity |

To rebuild the DIAMOND database if needed:
```bash
diamond makedb \
  --in data/uniref90.fasta \
  --db diamond_db/uniref90_diamond \
  --threads 32
```

### Sensitive search — jackhmmer

Uses [HMMER](http://hmmer.org/) `jackhmmer` for iterative profile-profile search
against `data/uniref90.fasta`. More sensitive than MMseqs2 for distant homologs
(e.g. <30% identity). Reads the FASTA sequentially, so it performs well on NFS.
Typical runtime: **3–15 minutes per query**.

- Installed at: `/usr/bin/jackhmmer` (HMMER 3.3.2)
- Does **not** require the MMseqs2 server to be running
- `identity` and `coverage` columns are `NaN` in results (not reported by jackhmmer)

---

## Python Query Interface

`uniprotdb.py` provides a high-level API. All methods return pandas DataFrames or
plain Python dicts.

```python
from uniprotdb import UniProtDB

db = UniProtDB(
    db_path="/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/uniprot.sqlite",
    diamond_db_path="/auto/sahara/namib/home/jhilzinger/databases/uniprot_db/diamond_db/uniref90_diamond.dmnd",
)

# Fast sequence search (DIAMOND) — 2–8 min depending on sensitivity
hits = db.search_by_sequence(fasta_str, mode="fast", sensitivity="sensitive", max_hits=100)

# Sensitive sequence search (jackhmmer) — 3–15 min, better for distant homologs
hits = db.search_by_sequence(fasta_str, mode="sensitive", jackhmmer_iterations=3)

# Look up a specific protein
info = db.get_by_accession("P04637")        # sequence, domains, features, cluster

# Find all proteins with a given domain
proteins = db.search_by_domain(pfam_id="PF00069")

# Get all members of a UniRef90 cluster
members = db.get_cluster_members("UniRef90_P04637")
```

Both modes return a DataFrame with identical columns:
`uniref90_id | rep_accession | evalue | identity | coverage | interpro_ids | pfam_ids | member_count | taxonomy_id`

(`identity` and `coverage` are `NaN` for jackhmmer results.)

---

## Known Issues

### MMseqs2 → DIAMOND migration (2026-03-03)

The original MMseqs2 index (`mmseqs_db/`, 619 GB) was removed due to incompatibility
with NFS-mounted storage (memory-mapped index files cause severe performance degradation
on NFS). DIAMOND replaces it as the fast search engine. DIAMOND performs sequential
disk reads and is fully NFS-compatible. The MMseqs2 index is not recoverable from
existing data but can be rebuilt with `mmseqs createdb` + `mmseqs createindex` if
local (non-NFS) storage becomes available.

---

### Domains table — re-parse (resolved 2026-03-09)

The original build assumed `protein2ipr.dat` had 7 columns but the actual file has 6
(no separate `db_name` column; the database name is encoded in the accession prefix).
This caused all rows to be rejected as malformed, leaving the `domains` table empty.

The bug in `parse_protein2ipr.py` was fixed on 2026-03-03. The re-parse completed
2026-03-09, inserting **1,174,000,000 domain rows**. The `domains` table is now fully
populated and all domain-related API calls are functional.

---

## Rebuilding / Updating

The build is fully checkpointed. Each step writes a `.step_complete` flag file and
is skipped on subsequent runs. To re-run a specific step, delete its flag:

```bash
# Re-run only the protein2ipr parsing step
rm /auto/sahara/namib/home/jhilzinger/databases/uniprot_db/.parse_protein2ipr_complete
bash build_database.sh
```

To rebuild everything from scratch, delete all flag files and the SQLite database.
Downloaded source files in `data/` are preserved and re-verified by checksum.
