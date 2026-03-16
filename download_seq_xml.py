"""Download patent sequence listings from USPTO bulk XML data.

For each patent, downloads the weekly APPXML zip containing the pgpub XML,
extracts the sequence-cwu block, parses sequences, and writes FASTA files.
Caches raw XML blocks locally to avoid re-downloading.
"""

import html
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
import tempfile
import time
import zipfile
from functools import partial
from pathlib import Path

# Unbuffered output
print = partial(print, flush=True)

import httpx
from dotenv import load_dotenv

from download_patents import read_patents, app_number_to_serial, safe_filename

load_dotenv()

DATA_DIR = Path("data/patents")
ASSOC_DOCS_PATH = Path("data/associated_docs.json")
PROGRESS_PATH = DATA_DIR / "seq_xml_progress.json"

API_BASE = "https://api.uspto.gov/api/v1"

# Three-letter to one-letter amino acid mapping
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Asx": "B", "Glx": "Z", "Xaa": "X",
    "Xle": "J",
}



# ---------------------------------------------------------------------------
# Progress tracking
# ---------------------------------------------------------------------------


def load_progress() -> dict:
    if PROGRESS_PATH.exists():
        return json.loads(PROGRESS_PATH.read_text())
    return {"done": [], "no_seq": [], "errors": []}


def save_progress(progress: dict):
    PROGRESS_PATH.write_text(json.dumps(progress, indent=2))


# ---------------------------------------------------------------------------
# Associated documents lookup
# ---------------------------------------------------------------------------


def load_associated_docs() -> dict:
    """Load cached associated-documents data, keyed by display_key."""
    if not ASSOC_DOCS_PATH.exists():
        raise FileNotFoundError(
            f"{ASSOC_DOCS_PATH} not found. Run query_assoc_docs first."
        )
    with open(ASSOC_DOCS_PATH) as f:
        docs = json.load(f)
    return {d["display_key"]: d for d in docs}


# ---------------------------------------------------------------------------
# Zip download and sequence extraction
# ---------------------------------------------------------------------------


def download_and_extract_sequences(
    doc: dict, client: httpx.Client
) -> list[dict] | None:
    """Download the APPXML zip, extract sequence-cwu block, parse sequences."""
    zip_name = doc["pgpub_zip"]
    xml_name = doc["pgpub_xml"]
    pub_number = xml_name.replace(".xml", "").split("_")[1]

    if not zip_name:
        return None

    # Extract year from zip name: ipa161013.zip -> 2016
    m = re.match(r"ipa(\d{2})(\d{4})", zip_name)
    if not m:
        print(f"    cannot parse zip name: {zip_name}")
        return None
    year = "20" + m.group(1)

    # Download zip via API (follows redirect to signed CloudFront URL)
    api_key = os.environ["USPTO_ODP_API_KEY"]
    url = f"{API_BASE}/datasets/products/files/APPXML/{year}/{zip_name}"

    print(f"    downloading {zip_name}...")
    resp = client.get(
        url,
        headers={"X-API-Key": api_key},
        follow_redirects=True,
        timeout=300,
    )

    if resp.status_code != 200 or len(resp.content) < 1000:
        print(f"    download failed: status {resp.status_code}")
        return None

    # Write to temp file and process
    with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp:
        tmp.write(resp.content)
        tmp_path = tmp.name

    try:
        return extract_sequences_from_zip(tmp_path, pub_number)
    finally:
        os.unlink(tmp_path)


def cached_cwu_path(patent_id: str) -> Path:
    """Return the path for a cached sequence-cwu XML file."""
    return DATA_DIR / f"{patent_id}_seq_cwu.xml"


def load_cached_cwu(patent_id: str) -> str | None:
    """Load a cached sequence-cwu XML block if it exists."""
    path = cached_cwu_path(patent_id)
    if path.exists():
        return path.read_text()
    return None


def extract_sequences_from_zip(
    zip_path: str, pub_number: str, patent_id: str | None = None,
) -> list[dict] | None:
    """Extract sequence-cwu block for a specific patent from an APPXML zip.

    If patent_id is provided, caches the raw XML block to disk.
    """
    with zipfile.ZipFile(zip_path) as z:
        # The zip contains one consolidated XML file
        xml_files = z.namelist()
        if not xml_files:
            return None

        content = z.read(xml_files[0]).decode("utf-8", errors="replace")

    # Find the sequence-cwu block for this patent's pub number
    # Each block: <sequence-cwu ...>...<doc-number>NNNNN</doc-number>...</sequence-cwu>
    starts = [m.start() for m in re.finditer(r"<sequence-cwu", content)]
    ends = [m.end() for m in re.finditer(r"</sequence-cwu>", content)]

    for s, e in zip(starts, ends):
        block = content[s:e]
        m = re.search(r"<doc-number>(\d+)</doc-number>", block[:500])
        if m and m.group(1) == pub_number:
            # Cache the raw XML block
            if patent_id:
                cached_cwu_path(patent_id).write_text(block)
            return parse_sequence_cwu(block)

    print(f"    pub {pub_number} not found in zip ({len(starts)} seq blocks)")
    return None


def parse_sequence_cwu(block: str) -> list[dict]:
    """Parse a sequence-cwu block into sequence records.

    Handles three formats:
    1. Modern <s210>/<s211>/<s212>/<s400> XML tags
    2. ST.26 XML with <INSDSeq_sequence> tags
    3. Legacy table-formatted ST.25 with HTML-encoded angle brackets
    """
    # Format 1: Modern XML tags
    if "<s210>" in block:
        return parse_s_tags(block)

    # Format 2: ST.26 XML
    if "<INSDSeq_sequence>" in block:
        return parse_st26(block)

    # Format 3: Legacy table format with HTML-encoded ST.25
    entries = re.findall(r"<entry>(.*?)</entry>", block)
    if entries:
        text = "\n".join(html.unescape(e) for e in entries)
        return parse_st25(text)

    return []


def parse_st26(block: str) -> list[dict]:
    """Parse ST.26 XML format with <SequenceData>/<INSDSeq> elements."""
    root = ET.fromstring(block)
    sequences = []

    for seq_data in root.iter("SequenceData"):
        seq_id_str = seq_data.get("sequenceIDNumber")
        if not seq_id_str:
            continue
        seq_id = int(seq_id_str)

        insd = seq_data.find("INSDSeq")
        if insd is None:
            continue

        len_el = insd.find("INSDSeq_length")
        type_el = insd.find("INSDSeq_moltype")
        seq_el = insd.find("INSDSeq_sequence")

        if seq_el is None or not seq_el.text:
            continue

        length = int(len_el.text) if len_el is not None and len_el.text else 0
        if length == 0:
            continue

        moltype = type_el.text.strip() if type_el is not None and type_el.text else ""
        seq_type = "PRT" if moltype == "AA" else "DNA"
        sequence = seq_el.text.strip()

        if sequence:
            sequences.append({
                "seq_id": seq_id,
                "length": length,
                "type": seq_type,
                "sequence": sequence,
            })

    return sequences


def parse_s_tags(block: str) -> list[dict]:
    """Parse modern <s210>/<s400> XML sequence format using an XML parser.

    Structure: <s200> contains header tags (s210 id, s211 length, s212 type),
    and sibling <s400> elements contain the sequence data.
    """
    root = ET.fromstring(block)
    sequences = []

    s200s = root.findall(".//s200")
    s400s = root.findall(".//s400")

    # Build s400 lookup keyed by sequence number
    s400_by_id = {}
    for s400 in s400s:
        text = s400.text or ""
        m = re.match(r"\s*(\d+)\s*\n", text)
        if m:
            s400_by_id[int(m.group(1))] = text[m.end():]

    for s200 in s200s:
        id_el = s200.find("s210")
        len_el = s200.find("s211")
        type_el = s200.find("s212")

        if id_el is None or not id_el.text:
            continue

        seq_id = int(id_el.text)
        length = int(len_el.text) if len_el is not None and len_el.text else 0
        raw_type = type_el.text.strip() if type_el is not None and type_el.text else ""

        if length == 0:
            continue

        raw_seq = s400_by_id.get(seq_id, "").strip()
        if not raw_seq:
            continue

        seq_type = "PRT" if raw_type == "PRT" else "DNA"

        if seq_type == "PRT":
            sequence = parse_protein_sequence(raw_seq)
        else:
            sequence = parse_dna_sequence(raw_seq)

        if sequence:
            sequences.append({
                "seq_id": seq_id,
                "length": length,
                "type": seq_type,
                "sequence": sequence,
            })

    return sequences


# ---------------------------------------------------------------------------
# ST.25 parsing (for legacy XML format with HTML-encoded entries)
# ---------------------------------------------------------------------------


def normalize_st25_tags(text: str) -> str:
    """Fix common artifacts in ST.25 angle-bracket tags."""
    text = re.sub(r"<\d*[Zz](\d{2,3})>", r"<\1>", text)
    text = re.sub(r"<[Zz]?\d?1[l1][l1][l1]?>", "<211>", text)
    text = re.sub(r"<[Zz]?\d?1[l1]2>", "<212>", text)
    text = re.sub(r"<[Zz]?\d?1[l1]3>", "<213>", text)
    text = re.sub(r"<[Zz]2[12]0>", "<220>", text)
    text = re.sub(r"<[Zz]2[12]1>", "<221>", text)
    text = re.sub(r"<[Zz]2[12]3>", "<223>", text)
    text = re.sub(r"<10>\s*SE[QOC] ID", "<210> SEQ ID", text)
    text = re.sub(r"<11>\s*LENGTH", "<211> LENGTH", text)
    text = re.sub(r"<12>\s*TYPE", "<212> TYPE", text)
    text = re.sub(r"<13>\s*ORGANISM", "<213> ORGANISM", text)
    text = re.sub(r"<210>\s*SE[QOC] ID NO[.:]*\s*", "<210> ", text)
    text = re.sub(r"<400>\s*SEQU?E?N?C?E?[.:]*\s*", "<400> ", text)
    text = re.sub(r"<211>\s*LENGTH[.:]*\s*", "<211> ", text)
    text = re.sub(r"<212>\s*TYPE[.:]*\s*", "<212> ", text)
    text = re.sub(r"<213>\s*ORGANISM[.:]*\s*", "<213> ", text)
    text = re.sub(r"(<[1-4]\d0>)\s*(\d+)/(\d+)", lambda m: m.group(1) + " " + m.group(2) + m.group(3), text)
    return text


def _infer_seq_type(raw_seq: str) -> str:
    """Infer sequence type from raw content."""
    if re.search(r"(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)\s", raw_seq):
        return "PRT"
    cleaned = re.sub(r"[\s\d]", "", raw_seq)
    if cleaned and all(c in "acgtun" for c in cleaned.lower()[:100]):
        return "DNA"
    return "UNK"


def parse_st25(text: str) -> list[dict]:
    """Parse ST.25 formatted text into a list of sequence records."""
    text = normalize_st25_tags(text)
    sequences = []
    blocks = re.split(r"<210>\s*", text)
    for block in blocks[1:]:
        record = {}
        m = re.match(r"(\d+)", block)
        if not m:
            continue
        record["seq_id"] = int(m.group(1))
        m = re.search(r"<211>\s*(\d+)", block)
        declared_length = int(m.group(1)) if m else None
        m = re.search(r"<212>\s*(DNA|RNA|PRT)", block)
        declared_type = m.group(1) if m else None
        if declared_length == 0 and m:
            continue
        m = re.search(r"<400>\s*\d+\s*\n(.*?)(?=<210>|\Z)", block, re.DOTALL)
        if not m:
            continue
        raw_seq = m.group(1).strip()
        if not raw_seq:
            continue
        seq_type = declared_type or _infer_seq_type(raw_seq)
        if seq_type == "PRT":
            record["sequence"] = parse_protein_sequence(raw_seq)
        elif seq_type in ("DNA", "RNA"):
            record["sequence"] = parse_dna_sequence(raw_seq)
        else:
            record["sequence"] = parse_dna_sequence(raw_seq)
            if not record["sequence"]:
                record["sequence"] = parse_protein_sequence(raw_seq)
            seq_type = "DNA" if record.get("sequence") and all(
                c in "acgtun" for c in record["sequence"][:50]
            ) else "PRT"
        record["type"] = seq_type
        record["length"] = declared_length or len(record.get("sequence", ""))
        if record.get("sequence"):
            sequences.append(record)
    return sequences


def parse_protein_sequence(raw: str) -> str:
    """Convert three-letter amino acid codes to single-letter FASTA sequence."""
    lines = raw.split("\n")
    cleaned_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if re.match(r"^[\d\s]+$", stripped):
            continue
        cleaned_lines.append(stripped)
    text = " ".join(cleaned_lines)
    one_letter = []
    tokens = re.findall(r"[A-Z][a-z]{1,3}", text)
    for token in tokens:
        if token in AA3_TO_1:
            one_letter.append(AA3_TO_1[token])
    return "".join(one_letter)


def parse_dna_sequence(raw: str) -> str:
    """Extract DNA/RNA sequence from ST.25 format."""
    seq = ""
    for line in raw.split("\n"):
        cleaned = re.sub(r"\d+\s*$", "", line)
        cleaned = re.sub(r"[^acgtunACGTUN]", "", cleaned)
        seq += cleaned.lower()
    return seq


def write_fasta_files(sequences: list[dict], patent_id: str, patent_number: str,
                      append: bool = False):
    """Write nucleotide and protein FASTA files."""
    nuc_seqs = []
    prot_seqs = []
    for seq in sequences:
        sid = seq["seq_id"]
        header = f">SEQ_{sid} Sequence {sid} from patent {patent_number}"
        entry = f"{header}\n{seq['sequence']}"
        if seq["type"] == "PRT":
            prot_seqs.append(entry)
        else:
            nuc_seqs.append(entry)
    nuc_path = DATA_DIR / f"{patent_id}_ocr_nuc.fasta"
    prot_path = DATA_DIR / f"{patent_id}_ocr_prot.fasta"
    mode = "a" if append else "w"
    if nuc_seqs:
        with open(nuc_path, mode) as fh:
            fh.write("\n\n".join(nuc_seqs) + "\n")
    if prot_seqs:
        with open(prot_path, mode) as fh:
            fh.write("\n\n".join(prot_seqs) + "\n")
    return len(nuc_seqs), len(prot_seqs)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    patents = read_patents()
    patents = [p for p in patents if p["seq_count"] > 0]
    assoc_docs = load_associated_docs()
    progress = load_progress()

    print(f"Processing {len(patents)} patents with sequences "
          f"({len(progress['done'])} already done)")

    # First pass: reparse any patents that have cached XML
    remaining = []
    done_or_skip = set(progress["done"] + progress["no_seq"])
    for patent in patents:
        dk = patent["display_key"]
        pid = safe_filename(dk)
        if pid in done_or_skip:
            continue

        m_pn = re.match(r"(US\s+\S+)", dk)
        patent_number = m_pn.group(1) if m_pn else dk

        cached = load_cached_cwu(pid)
        if cached:
            print(f"  {dk} (cached)...", end=" ")
            sequences = parse_sequence_cwu(cached)
            if sequences:
                n_nuc, n_prot = write_fasta_files(sequences, pid, patent_number)
                total = n_nuc + n_prot
                print(f"{total} sequences ({n_nuc} nuc + {n_prot} prot)")
                progress["done"].append(pid)
            else:
                print("no sequences found")
                if pid not in progress["no_seq"]:
                    progress["no_seq"].append(pid)
            save_progress(progress)
        else:
            remaining.append(patent)

    # Group remaining patents by zip to avoid downloading the same zip twice
    by_zip = {}
    done_or_skip = set(progress["done"] + progress["no_seq"])
    for patent in remaining:
        dk = patent["display_key"]
        pid = safe_filename(dk)
        if pid in done_or_skip:
            continue
        doc = assoc_docs.get(dk)
        if not doc:
            continue
        zip_name = doc.get("pgpub_zip", "")
        if not zip_name:
            continue
        by_zip.setdefault(zip_name, []).append((patent, doc))

    if by_zip:
        print(f"\nNeed to download {len(by_zip)} unique zips for "
              f"{sum(len(v) for v in by_zip.values())} patents")

    with httpx.Client() as client:
        for i, (zip_name, patent_docs) in enumerate(sorted(by_zip.items())):
            print(f"\n[{i + 1}/{len(by_zip)}] {zip_name} "
                  f"({len(patent_docs)} patents)")

            # Download zip once
            api_key = os.environ["USPTO_ODP_API_KEY"]
            m = re.match(r"ipa(\d{2})(\d{4})", zip_name)
            if not m:
                continue
            year = "20" + m.group(1)
            url = f"{API_BASE}/datasets/products/files/APPXML/{year}/{zip_name}"

            print(f"  downloading {zip_name}...")
            try:
                resp = client.get(
                    url,
                    headers={"X-API-Key": api_key},
                    follow_redirects=True,
                    timeout=300,
                )
            except Exception as e:
                print(f"  download error: {e}")
                for patent, doc in patent_docs:
                    pid = safe_filename(patent["display_key"])
                    if pid not in progress["errors"]:
                        progress["errors"].append(pid)
                save_progress(progress)
                continue

            if resp.status_code != 200 or len(resp.content) < 1000:
                print(f"  download failed: status {resp.status_code}")
                continue

            zip_size_mb = len(resp.content) / (1024 * 1024)
            print(f"  downloaded {zip_size_mb:.1f} MB")

            # Write to temp file
            with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp:
                tmp.write(resp.content)
                tmp_path = tmp.name

            try:
                for patent, doc in patent_docs:
                    dk = patent["display_key"]
                    pid = safe_filename(dk)
                    pub_number = doc["pgpub_xml"].replace(".xml", "").split("_")[1]

                    m_pn = re.match(r"(US\s+\S+)", dk)
                    patent_number = m_pn.group(1) if m_pn else dk

                    print(f"  {dk}...", end=" ")

                    sequences = extract_sequences_from_zip(
                        tmp_path, pub_number, patent_id=pid,
                    )

                    if sequences:
                        n_nuc, n_prot = write_fasta_files(
                            sequences, pid, patent_number
                        )
                        total = n_nuc + n_prot
                        print(f"{total} sequences ({n_nuc} nuc + {n_prot} prot)")
                        progress["done"].append(pid)
                    else:
                        print("no sequences found")
                        if pid not in progress["no_seq"]:
                            progress["no_seq"].append(pid)

                    save_progress(progress)
            finally:
                os.unlink(tmp_path)

            time.sleep(1)  # Be polite to the API

    print(f"\n{'=' * 50}")
    print(f"Done:    {len(progress['done'])}/{len(patents)}")
    print(f"No seq:  {len(progress['no_seq'])}")
    print(f"Errors:  {len(progress['errors'])}")


if __name__ == "__main__":
    main()
