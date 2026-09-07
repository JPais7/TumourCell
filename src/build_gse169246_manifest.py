#!/usr/bin/env python3
"""Build a sample manifest for GSE169246 without reading expression counts."""

from __future__ import annotations

import csv
import gzip
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOFT = ROOT / "tmp/GSE169246/GSE169246_family.soft.gz"
BARCODES = ROOT / "data/raw/GSE169246/GSE169246_TNBC_RNA.barcode.tsv.gz"
OUTPUT = ROOT / "data/manifests/GSE169246_RNA_sample_manifest.csv"
GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169246"
BARCODE_SHA256 = "c61d02185e54905a1f7ed658d69fd588db2526a9130951272240338c083a9141"


def parse_soft() -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(SOFT, "rt", encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    records.append(current)
                current = {"geo_accession": line.split(" = ", 1)[1]}
            elif current is not None and line.startswith("!Sample_title = "):
                current["sample_title"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_source_name_ch1 = "):
                current["patient_id"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_library_strategy = "):
                current["library_strategy"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key.replace(" ", "_")] = item
    if current:
        records.append(current)
    return records


def barcode_counts() -> Counter[str]:
    counts: Counter[str] = Counter()
    with gzip.open(BARCODES, "rt", encoding="utf-8") as stream:
        for line in stream:
            counts[line.rstrip("\n").rsplit(".", 1)[1]] += 1
    return counts


def normalize_regimen(value: str) -> str:
    compact = value.lower().replace("-", "")
    if compact == "chemo":
        return "paclitaxel"
    if "pdl1" in compact and "chemo" in compact:
        return "paclitaxel+atezolizumab"
    raise ValueError(f"unexpected treatment label: {value}")


def main() -> None:
    counts = barcode_counts()
    rows = []
    for record in parse_soft():
        title = record.get("sample_title", "")
        if record.get("library_strategy") != "RNA-Seq" or title not in counts:
            continue
        timepoint, patient, compartment = title.split("_")
        rows.append(
            {
                "study_id": "GSE169246",
                "geo_accession": record["geo_accession"],
                "sample_title": title,
                "patient_id": patient,
                "timepoint": timepoint,
                "compartment": {"b": "blood", "t": "tumour"}[compartment],
                "anatomical_site": record.get("tissue", ""),
                "regimen": normalize_regimen(record.get("treatment", "")),
                "regimen_source": record.get("treatment", ""),
                "response_primary_source": "",
                "cd45_enriched": "true",
                "cell_barcodes": counts[title],
                "source_url": GEO_URL,
                "barcode_file_sha256": BARCODE_SHA256,
            }
        )

    if len(rows) != len(counts):
        missing = sorted(set(counts) - {row["sample_title"] for row in rows})
        raise ValueError(f"barcode samples absent from RNA-Seq SOFT records: {missing}")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} RNA samples and {sum(counts.values())} barcodes -> {OUTPUT}")


if __name__ == "__main__":
    main()
