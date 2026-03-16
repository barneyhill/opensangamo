"""Download patent PDFs for all patents in lens-export.csv.

PDFs: image-ppubs.uspto.gov (no auth required)
"""

import csv
import json
import re
import time
from pathlib import Path

import httpx

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

PDF_URL = "https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/{number}"

DATA_DIR = Path("data/patents")
CSV_PATH = Path("data/lens-export.csv")
PROGRESS_PATH = DATA_DIR / "progress.json"

# ---------------------------------------------------------------------------
# Progress tracking
# ---------------------------------------------------------------------------


def load_progress() -> dict:
    if PROGRESS_PATH.exists():
        return json.loads(PROGRESS_PATH.read_text())
    return {"pdf_done": [], "seq_done": [], "seq_not_found": []}


def save_progress(progress: dict):
    PROGRESS_PATH.write_text(json.dumps(progress, indent=2))


# ---------------------------------------------------------------------------
# Patent number parsing
# ---------------------------------------------------------------------------


def display_key_to_number(display_key: str) -> str:
    """Extract the bare number used by image-ppubs.

    'US 9937207 B2'      -> '9937207'
    'US 2019/0046579 A1' -> '20190046579'
    """
    key = display_key.strip()
    m = re.match(r"US\s+(\d{4})/(\d{7})\s+A\d", key)
    if m:
        return m.group(1) + m.group(2)
    m = re.match(r"US\s+(\d+)\s+B\d", key)
    if m:
        return m.group(1)
    raise ValueError(f"Cannot parse display key: {display_key}")


def app_number_to_serial(app_number: str) -> str | None:
    """Extract application serial from CSV format.

    'US 201414221074 A' -> '14221074'
    """
    m = re.match(r"US\s+\d{4}(\d+)\s+A", app_number.strip())
    return m.group(1) if m else None


def safe_filename(display_key: str) -> str:
    return display_key.replace(" ", "_").replace("/", "-")


# ---------------------------------------------------------------------------
# CSV reading
# ---------------------------------------------------------------------------


def read_patents() -> list[dict]:
    """Read patent list from CSV, sorted by sequence count descending."""
    patents = []
    with open(CSV_PATH, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            patents.append({
                "display_key": row["Display Key"],
                "app_number": row["Application Number"],
                "title": row["Title"],
                "seq_count": int(row["Sequence Count"] or 0),
            })
    patents.sort(key=lambda p: p["seq_count"], reverse=True)
    return patents


# ---------------------------------------------------------------------------
# PDF download
# ---------------------------------------------------------------------------


def download_pdf(patent: dict, client: httpx.Client) -> bool:
    """Download full patent PDF from image-ppubs (no auth)."""
    number = display_key_to_number(patent["display_key"])
    path = DATA_DIR / f"{safe_filename(patent['display_key'])}.pdf"
    if path.exists():
        return True

    url = PDF_URL.format(number=number)
    print(f"  PDF: {url}")
    resp = client.get(url, follow_redirects=True, timeout=60)
    if resp.status_code == 200 and len(resp.content) > 1000:
        path.write_bytes(resp.content)
        print(f"  PDF: saved {path.name} ({len(resp.content) // 1024} KB)")
        return True

    print(f"  PDF: failed (status={resp.status_code})")
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    patents = read_patents()
    progress = load_progress()
    print(f"Loaded {len(patents)} patents ({len(progress['pdf_done'])} PDFs cached)")

    with httpx.Client() as client:
        for i, patent in enumerate(patents):
            dk = patent["display_key"]
            print(f"\n[{i + 1}/{len(patents)}] {dk}"
                  f" — {patent['title'][:55]}... ({patent['seq_count']} seqs)")

            try:
                if dk not in progress["pdf_done"]:
                    if download_pdf(patent, client):
                        progress["pdf_done"].append(dk)
                        save_progress(progress)
            except Exception as e:
                print(f"  ERROR: {e}")
                continue

    n = len(patents)
    print(f"\n{'=' * 40}")
    print(f"PDFs:      {len(progress['pdf_done'])}/{n}")


if __name__ == "__main__":
    main()
