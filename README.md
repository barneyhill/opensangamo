# opensangamo

A dataset of Sangamo Therapeutics zinc finger protein (ZFP) designs extracted from USPTO patents.

The final dataset is **`data/sangamo_zinc_fingers.csv`** — 1094 unique ZFP designs with target DNA sequences and recognition helix amino acid sequences (F1-F6).

![Screenshot](Screenshot.png)

## Dataset columns

| Column | Description |
|--------|-------------|
| `patent_id` | USPTO patent identifier |
| `ZFN` | Zinc finger nuclease / ZFP-TF identifier |
| `Target Sequence` | DNA target sequence (nucleotide) |
| `F1`-`F6` | Recognition helix amino acid sequences for fingers 1-6 |

## How the dataset was built

### 1. Patent discovery

`data/lens-export.csv` — 115 Sangamo patents with sequence listings, exported from [Lens.org](https://www.lens.org/lens/search/patent/list?applicant.must=SANGAMO%20THERAPEUTICS%20INC).

### 2. Download patent PDFs

```
uv run python download_patents.py
```

Downloads full patent PDFs from `image-ppubs.uspto.gov` into `data/patents/`.

### 3. Extract ZFP tables via LLM

```
uv run python server.py
```

FastAPI web app for browsing patent PDFs, cropping table regions, and running LLM-powered extraction. Outputs `_extractions.json` files containing ZFN IDs, target sequences, and F1-F6 helix references (as SEQ ID numbers).

### 4. Download sequence listings from USPTO bulk XML

```
uv run python download_seq_xml.py
```

For each patent, downloads the weekly APPXML zip from `api.uspto.gov`, extracts the `<sequence-cwu>` block, and parses sequences (ST.25 and ST.26 formats). Writes FASTA files (`_ocr_nuc.fasta`, `_ocr_prot.fasta`) and caches raw XML blocks (`_seq_cwu.xml`) to avoid re-downloading.

### 5. Combine into final CSV

```
uv run python combine_extractions.py
```

Resolves SEQ ID references in the extraction JSONs to actual sequences using the FASTA files, then:
- Removes rows missing a target sequence or all finger columns
- Deduplicates on `(patent_id, ZFN)`
- Deduplicates on `(Target Sequence, F1..F6)` across patents

## Requirements

```
uv sync
```

Requires a `USPTO_ODP_API_KEY` in `.env` for sequence downloads (step 4).
