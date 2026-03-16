"""Combine patent extraction JSONs into a unified CSV dataset."""

import csv
import json
import glob
import re
import sys
from collections import defaultdict
from pathlib import Path

DATA_DIR = Path("data/patents")
OUTPUT_CSV = Path("data/sangamo_zinc_fingers.csv")

# Column name mappings
ID_COLUMNS = {
    "ZFP-TF Name", "ZFP", "ZFP ID", "ZFN ID", "ZFN Name", "SBS #", "SBS ID",
    "ZFP #", "SBS", "SBS#", "ZFN", "ZFN Number", "ZFN OD", "SBS Number",
}
TARGET_COLUMNS = {
    "Binding Sequence", "Target Sequence", "Target Site", "Target site",
    "Target ID", "Target sequence",
}
F_COLUMNS = {"F1", "F2", "F3", "F4", "F5", "F6"}
AA_COLUMNS = {"AA Sequence"}
IGNORE_COLUMNS = set()

OUTPUT_FIELDS = ["patent_id", "ZFN", "Target Sequence", "F1", "F2", "F3", "F4", "F5", "F6", "AA Sequence"]

warnings = []


def warn(msg):
    warnings.append(msg)


def parse_fasta_files(patent_id):
    """Parse nuc and prot FASTA files, return {seq_id_int: sequence}."""
    seqs = {}
    suffixes = [
        "_ocr_nuc.fasta", "_ocr_prot.fasta",
    ]
    found_any = False
    for suffix in suffixes:
        path = DATA_DIR / f"{patent_id}{suffix}"
        if not path.exists():
            continue
        found_any = True
        with open(path) as fh:
            content = fh.read()
        for m in re.finditer(
            r">(\S+)\s+Sequence (\d+) from patent[^\n]*\n([^>]+)", content
        ):
            seq_id = int(m.group(2))
            seq = m.group(3).replace("\n", "").strip()
            seqs[seq_id] = seq
    if not found_any:
        warn(f"No FASTA files found for {patent_id}")
    return seqs


def extract_zfn_id(raw):
    """Extract integer ID from heterogeneous ZFN ID strings."""
    if raw is None:
        return None
    raw = str(raw).strip()
    if not raw:
        return None
    # Try to extract digit sequence
    digits = re.findall(r"\d+", raw)
    if digits:
        # Use the longest digit sequence
        return max(digits, key=len)
    # No digits found — keep original, warn
    warn(f"Non-numeric ZFN ID kept as-is: {raw!r}")
    return raw


def resolve_value(val, seq_lookup, field_name, patent_id):
    """Resolve a SEQ ID integer (or other value) to its actual sequence."""
    if val is None or val == "" or val == 0:
        return ""
    # Handle lists — flatten to first element
    if isinstance(val, list):
        if len(val) == 0:
            return ""
        if len(val) == 1:
            val = val[0]
        else:
            # Multiple values — use first, warn
            warn(f"{patent_id} {field_name}: array {val}, using first element")
            val = val[0]
    # Handle string values
    if isinstance(val, str):
        val = val.strip()
        if val in ("N/A", "n/a", "NA", "-", ""):
            return ""
        # Try to extract SEQ ID from strings like "QKINLVN (SEQ ID NO: 61)"
        m = re.search(r"SEQ ID (?:NO[.:]?\s*)?(\d+)", val, re.IGNORECASE)
        if m:
            seq_id = int(m.group(1))
            if seq_id in seq_lookup:
                return seq_lookup[seq_id]
            warn(f"{patent_id} {field_name}: SEQ ID {seq_id} from string not in FASTA")
            return ""
        # If it's a pure number string, treat as seq ID
        if val.isdigit():
            val = int(val)
        else:
            # Already a sequence string or other text — keep as-is
            return val
    # Integer SEQ ID — look up
    if isinstance(val, (int, float)):
        seq_id = int(val)
        if seq_id == 0:
            return ""
        if seq_id in seq_lookup:
            return seq_lookup[seq_id]
        warn(f"{patent_id} {field_name}: SEQ ID {seq_id} not found in FASTA")
        return ""
    return str(val)


_NUC_RE = re.compile(r"^[ACGTacgtNn]+$")
_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYO]+$")


def validate_sequence(val, field_name, patent_id):
    """Blank out values that aren't the expected type (nuc for target, aa for F*)."""
    if not val:
        return ""
    if field_name == "Target Sequence":
        if _NUC_RE.match(val):
            return val
        warn(f"{patent_id} {field_name}: not nucleotide, clearing: {val[:40]}")
        return ""
    # F1-F6
    if _AA_RE.match(val):
        return val
    warn(f"{patent_id} {field_name}: not amino acid, clearing: {val[:40]}")
    return ""


def process_patent(patent_id, extractions, seq_lookup):
    """Process all extraction tables for a patent, return list of normalized rows."""
    rows_by_zfn = defaultdict(dict)  # zfn_id -> merged fields

    for table in extractions:
        columns = table.get("columns", [])
        col_headers = [c["header"] for c in columns]

        # Identify column roles
        id_col = None
        target_col = None
        aa_col = None
        f_cols = {}
        for header in col_headers:
            if header in ID_COLUMNS:
                id_col = header
            elif header in TARGET_COLUMNS:
                target_col = header
            elif header in F_COLUMNS:
                f_cols[header] = header
            elif header in AA_COLUMNS:
                aa_col = header
            elif header not in IGNORE_COLUMNS:
                warn(f"{patent_id}: unknown column {header!r}")

        if id_col is None:
            warn(f"{patent_id} page {table.get('page')}: no ID column found in {col_headers}")
            continue

        for row in table.get("rows", []):
            zfn = extract_zfn_id(row.get(id_col))
            if zfn is None:
                continue

            existing = rows_by_zfn[zfn]

            # Merge target
            if target_col and target_col in row:
                val = resolve_value(row[target_col], seq_lookup, "Target Sequence", patent_id)
                val = validate_sequence(val, "Target Sequence", patent_id)
                if val and not existing.get("Target Sequence"):
                    existing["Target Sequence"] = val

            # Merge F1-F6
            for f in ["F1", "F2", "F3", "F4", "F5", "F6"]:
                if f in row:
                    val = resolve_value(row[f], seq_lookup, f, patent_id)
                    val = validate_sequence(val, f, patent_id)
                    if val and not existing.get(f):
                        existing[f] = val

            # Merge AA Sequence
            if aa_col and aa_col in row:
                val = resolve_value(row[aa_col], seq_lookup, "AA Sequence", patent_id)
                if val and _AA_RE.match(val) and not existing.get("AA Sequence"):
                    existing["AA Sequence"] = val

    # Build output rows
    result = []
    for zfn, fields in rows_by_zfn.items():
        out = {"patent_id": patent_id, "ZFN": zfn}
        for col in ["Target Sequence", "F1", "F2", "F3", "F4", "F5", "F6", "AA Sequence"]:
            out[col] = fields.get(col, "")
        result.append(out)
    return result


def main():
    # Step 1: Find all patents
    extraction_files = sorted(glob.glob(str(DATA_DIR / "*_extractions.json")))
    print(f"Found {len(extraction_files)} extraction files")

    # Step 2: Parse FASTA files
    patent_ids = set()
    for f in extraction_files:
        pid = Path(f).name.replace("_extractions.json", "")
        patent_ids.add(pid)

    seq_lookups = {}
    for pid in sorted(patent_ids):
        seq_lookups[pid] = parse_fasta_files(pid)

    # Step 3: Process each patent
    all_rows = []
    per_patent_counts = {}

    for f in extraction_files:
        pid = Path(f).name.replace("_extractions.json", "")
        with open(f) as fh:
            extractions = json.load(fh)
        rows = process_patent(pid, extractions, seq_lookups.get(pid, {}))
        all_rows.extend(rows)
        per_patent_counts[pid] = len(rows)

    print(f"Total rows before filtering: {len(all_rows)}")

    # Step 4a: Remove rows with missing target seq or no F columns
    filtered = []
    removed_incomplete = 0
    for row in all_rows:
        target = row.get("Target Sequence", "").strip()
        has_any_f = any(row.get(f, "").strip() for f in ["F1", "F2", "F3", "F4", "F5", "F6"])
        if target and has_any_f:
            filtered.append(row)
        else:
            removed_incomplete += 1
    print(f"Removed incomplete (no target or no F*): {removed_incomplete}")

    # Step 4b: Dedup on (patent_id, ZFN)
    seen = {}
    deduped_zfn = []
    dup_zfn = 0
    for row in filtered:
        key = (row["patent_id"], row["ZFN"])
        if key in seen:
            prev = seen[key]
            is_identical = all(row.get(c) == prev.get(c) for c in OUTPUT_FIELDS[2:])
            if not is_identical:
                warn(f"Conflicting duplicate: {key}, keeping first")
            dup_zfn += 1
            continue
        seen[key] = row
        deduped_zfn.append(row)
    print(f"Duplicates removed (patent_id, ZFN): {dup_zfn}")

    # Step 4c: Dedup on (Target Sequence, F1..F6)
    seen_seq = {}
    deduped = []
    dup_seq = 0
    for row in deduped_zfn:
        key = tuple(row.get(c, "") for c in ["Target Sequence", "F1", "F2", "F3", "F4", "F5", "F6"])
        if key in seen_seq:
            dup_seq += 1
            continue
        seen_seq[key] = row
        deduped.append(row)
    print(f"Duplicates removed (target+fingers): {dup_seq}")

    print(f"Final row count: {len(deduped)}")

    # Step 5: Write CSV
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(deduped)

    print(f"Written to {OUTPUT_CSV}")

    # Summary
    print(f"\nPer-patent counts:")
    for pid in sorted(per_patent_counts):
        print(f"  {pid}: {per_patent_counts[pid]}")

    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for w in warnings:
            print(f"  {w}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
