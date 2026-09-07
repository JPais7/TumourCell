#!/usr/bin/env python3
"""Export the already-frozen P8 definition for Phase 2, without external data."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/phase1/blind_discovery/program_definitions_v1.npz"
OUTPUT = ROOT / "results/phase2/p8_definition_frozen.json"
EXPECTED = "c353b2b5edb47c18e960aeb10622dc673e6f6fdf24342deebd784a3982371aec"


def main() -> None:
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    if digest != EXPECTED:
        raise RuntimeError(f"Phase 1 definition hash changed: {digest}")
    frozen = np.load(SOURCE)
    genes = frozen["genes"].astype(str)
    weights = frozen["weights"][7].astype(float)
    scales = frozen["gene_scale"].astype(float)
    payload = {
        "phase1_definition_sha256": digest,
        "program_id": "P8",
        "program_label": "antigen-presentation-interferon",
        "rank": int(frozen["rank"]),
        "seed": int(frozen["seed"]),
        "gene_count": len(genes),
        "genes": genes.tolist(),
        "weights": weights.tolist(),
        "gene_scale": scales.tolist(),
        "score": "sum(log1p(CPM_gene) / discovery_gene_sd * P8_weight); P8 weights sum to one",
        "missing_gene_rule": "set missing contribution to zero and renormalize P8 weights over present genes only if coverage thresholds pass",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"{len(genes)} frozen genes -> {OUTPUT}")


if __name__ == "__main__":
    main()
