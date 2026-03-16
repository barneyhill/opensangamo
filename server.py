"""Patent table extraction web app — FastAPI backend."""

import csv
import json
import logging
import re
from datetime import datetime, timezone
from pathlib import Path

from dotenv import load_dotenv
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from openai import AsyncOpenAI
from pydantic import BaseModel

load_dotenv()

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

DATA_DIR = Path("data/patents")
CSV_PATH = Path("data/lens-export.csv")

log = logging.getLogger("extract")
logging.basicConfig(level=logging.INFO)

app = FastAPI()
_openai_client: AsyncOpenAI | None = None


def get_openai() -> AsyncOpenAI:
    global _openai_client
    if _openai_client is None:
        _openai_client = AsyncOpenAI()
    return _openai_client

# ---------------------------------------------------------------------------
# Patent helpers
# ---------------------------------------------------------------------------


def safe_filename(display_key: str) -> str:
    return display_key.replace(" ", "_").replace("/", "-")


def read_patents() -> list[dict]:
    """Read patent list from CSV, sorted by sequence count descending."""
    patents = []
    with open(CSV_PATH, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            dk = row["Display Key"]
            fid = safe_filename(dk)
            pdf_path = DATA_DIR / f"{fid}.pdf"
            if not pdf_path.exists():
                continue
            seq_count = int(row["Sequence Count"] or 0)
            extractions_path = DATA_DIR / f"{fid}_extractions.json"
            patents.append({
                "id": fid,
                "display_key": dk,
                "title": row["Title"],
                "seq_count": seq_count,
                "has_extractions": extractions_path.exists(),
            })
    patents.sort(key=lambda p: p["seq_count"], reverse=True)
    return patents


_patents_cache: list[dict] | None = None


def get_patents() -> list[dict]:
    global _patents_cache
    if _patents_cache is None:
        _patents_cache = read_patents()
    return _patents_cache


def invalidate_cache():
    global _patents_cache
    _patents_cache = None


# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------


@app.get("/")
async def index():
    return FileResponse("static/index.html")


@app.get("/api/patents")
async def list_patents():
    return get_patents()


@app.get("/api/patents/{patent_id}/pdf")
async def get_pdf(patent_id: str):
    path = DATA_DIR / f"{patent_id}.pdf"
    if not path.exists():
        raise HTTPException(404, "PDF not found")
    return FileResponse(path, media_type="application/pdf")


@app.get("/api/patents/{patent_id}/extractions")
async def get_extractions(patent_id: str):
    path = DATA_DIR / f"{patent_id}_extractions.json"
    if not path.exists():
        return []
    return json.loads(path.read_text())


class ColumnConfig(BaseModel):
    header: str
    mode: str  # "all" or "seq_ids"


class ExtractRequest(BaseModel):
    image_base64: str
    columns: list[ColumnConfig]


@app.post("/api/patents/{patent_id}/extract")
async def extract_table(patent_id: str, req: ExtractRequest):
    if not req.columns:
        raise HTTPException(400, "At least one column required")

    # Build column descriptions for the prompt
    col_descriptions = []
    for col in req.columns:
        if col.mode == "seq_ids":
            col_descriptions.append(
                f'"{col.header}" (SEQ_ID mode): cells contain sequences with '
                f"a parenthesised SEQ ID number like 'ATCG... (213)'. "
                f"Extract ONLY the integer from the parentheses. "
                f"If multiple SEQ IDs, return as an array of integers. "
                f"Always return integers, never strings."
            )
        else:
            col_descriptions.append(
                f'"{col.header}" (text mode): extract the full cell text as a string.'
            )

    prompt = (
        "Extract data from this table image. Return a JSON object with a single "
        '"rows" key containing an array of row objects.\n\n'
        "Columns to extract:\n"
        + "\n".join(f"- {d}" for d in col_descriptions)
        + "\n\nReturn JSON: {\"rows\": [{\"col_header\": value, ...}]}"
    )

    try:
        response = await get_openai().chat.completions.create(
            model="gpt-5-mini",
            messages=[
                {
                    "role": "user",
                    "content": [
                        {"type": "text", "text": prompt},
                        {
                            "type": "image_url",
                            "image_url": {
                                "url": f"data:image/png;base64,{req.image_base64}"
                            },
                        },
                    ],
                }
            ],
            response_format={"type": "json_object"},
        )
    except Exception as e:
        raise HTTPException(502, f"OpenAI API error: {e}")

    raw = response.choices[0].message.content
    log.info("Model response (first 500 chars): %s", raw[:500])
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        log.error("Invalid JSON from model: %s", raw[:500])
        raise HTTPException(502, f"Invalid JSON from model: {raw[:200]}")

    rows = data.get("rows", [])
    log.info("Extracted %d rows for patent %s", len(rows), patent_id)

    # Coerce SEQ ID columns to int / list[int]
    seq_cols = {col.header for col in req.columns if col.mode == "seq_ids"}
    for row in rows:
        for col_name in seq_cols:
            val = row.get(col_name)
            if val is None:
                continue
            row[col_name] = _coerce_seq_id(val)

    return {"rows": rows}


def _coerce_seq_id(val):
    """Coerce a value to int or list[int] for SEQ ID columns."""
    if isinstance(val, int):
        return val
    if isinstance(val, float):
        return int(val)
    if isinstance(val, str):
        # Try to extract integers from parentheses or bare numbers
        nums = re.findall(r"\d+", val)
        if len(nums) == 1:
            return int(nums[0])
        if len(nums) > 1:
            return [int(n) for n in nums]
        return val
    if isinstance(val, list):
        return [_coerce_seq_id(v) for v in val]
    return val


class SaveRequest(BaseModel):
    page: int
    crop: dict
    columns: list[ColumnConfig]
    rows: list[dict]


@app.post("/api/patents/{patent_id}/extractions")
async def save_extraction(patent_id: str, req: SaveRequest):
    path = DATA_DIR / f"{patent_id}_extractions.json"
    existing = []
    if path.exists():
        existing = json.loads(path.read_text())

    entry = {
        "page": req.page,
        "crop": req.crop,
        "columns": [c.model_dump() for c in req.columns],
        "rows": req.rows,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    existing.append(entry)
    path.write_text(json.dumps(existing, indent=2))
    invalidate_cache()
    return {"ok": True, "count": len(existing)}


@app.delete("/api/patents/{patent_id}/extractions/{index}")
async def delete_extraction(patent_id: str, index: int):
    path = DATA_DIR / f"{patent_id}_extractions.json"
    if not path.exists():
        raise HTTPException(404, "No extractions found")
    existing = json.loads(path.read_text())
    if index < 0 or index >= len(existing):
        raise HTTPException(404, "Extraction index out of range")
    existing.pop(index)
    path.write_text(json.dumps(existing, indent=2))
    invalidate_cache()
    return {"ok": True, "count": len(existing)}


# Static files (after API routes so they don't shadow)
app.mount("/static", StaticFiles(directory="static"), name="static")
