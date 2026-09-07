#!/usr/bin/env python3
"""Create auditable sample manifests from GEO SOFT records.

This script downloads metadata only. It does not download expression matrices.
"""

from __future__ import annotations

import csv
import gzip
import io
import re
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "manifests"

SERIES = {
    "GSE205472": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205472/soft/GSE205472_family.soft.gz",
    "GSE246613": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246613/soft/GSE246613_family.soft.gz",
}


def read_soft(url: str) -> str:
    with urllib.request.urlopen(url) as response:
        return gzip.decompress(response.read()).decode("utf-8", errors="replace")


def samples_from_soft(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE = "):
            if current:
                rows.append(current)
            current = {"gsm": line.split(" = ", 1)[1]}
        elif current is not None and line.startswith("!Sample_title = "):
            current["title"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_source_name_ch1 = "):
            current["source_name"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_organism_ch1 = "):
            current["organism"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
            value = line.split(" = ", 1)[1]
            if ": " in value:
                key, val = value.split(": ", 1)
                current[key.strip().lower().replace(" ", "_")] = val.strip()
        elif current is not None and line.startswith("!Sample_relation = SRA: "):
            current["sra_url"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_relation = BioSample: "):
            current["biosample_url"] = line.split(" = ", 1)[1]
    if current:
        rows.append(current)
    return rows


def normalize(series: str, sample: dict[str, str]) -> dict[str, str]:
    title = sample.get("title", "")
    if series == "GSE205472":
        patient = sample.get("source_name", "").replace("SC-", "")
        assay = "scRNA-seq"
        timepoint = sample.get("treatment_state", "")
        biopsy = "pre" if timepoint.lower().startswith("pre") else "post"
    else:
        match = re.match(r"(h\d+)([A-C])([0-9]*?)_(P|TCR|N)$", title)
        patient = match.group(1) if match else ""
        letter = match.group(2) if match else ""
        replicate = match.group(3) if match else ""
        assay = "scRNA-seq" if title.endswith("_P") else ("snRNA-seq" if title.endswith("_N") else "scTCR-seq")
        treatment = sample.get("treatment", "")
        biopsy = {"A": "baseline", "B": "post_pembrolizumab", "C": "post_pembrolizumab_radiation"}.get(letter, treatment)
        timepoint = biopsy + (("_rep" + replicate) if replicate else "")
    return {
        "study_id": series,
        "organism": sample.get("organism", ""),
        "analysis_scope": "human_clinical" if sample.get("organism") == "Homo sapiens" else "preclinical",
        "patient_id": patient,
        "sample_id": sample.get("gsm", ""),
        "sample_title": title,
        "biopsy_timepoint": timepoint,
        "treatment_label": sample.get("treatment_state", sample.get("treatment", "")),
        "response": sample.get("response", ""),
        "subtype": sample.get("disease_state", ""),
        "tissue": sample.get("tissue", sample.get("source_name", "")),
        "assay": assay,
        "biosample_url": sample.get("biosample_url", ""),
        "sra_url": sample.get("sra_url", ""),
        "source_soft_url": SERIES[series],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    columns = [
        "study_id", "organism", "analysis_scope", "patient_id", "sample_id", "sample_title", "biopsy_timepoint",
        "treatment_label", "response", "subtype", "tissue", "assay",
        "biosample_url", "sra_url", "source_soft_url",
    ]
    for series, url in SERIES.items():
        records = [normalize(series, item) for item in samples_from_soft(read_soft(url))]
        path = OUT / f"{series}_sample_manifest.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(records)
        patients = {row["patient_id"] for row in records if row["patient_id"]}
        human = sum(row["analysis_scope"] == "human_clinical" for row in records)
        print(f"{series}: {len(records)} assay records ({human} human), {len(patients)} parsed patients -> {path}")


if __name__ == "__main__":
    main()
