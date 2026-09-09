"""Deterministic, versioned encoder for extensible continuous programme scores."""
from __future__ import annotations
from dataclasses import dataclass
import hashlib, json
from pathlib import Path
from typing import Mapping, Sequence
import numpy as np
import yaml
from .latent_state import LatentFeature, LatentState
from .uncertainty import uncertainty_components, quality_label


@dataclass(frozen=True)
class ProgramSpec:
    name: str; genes: tuple[str, ...]; weights: tuple[float, ...]; role: str = "candidate_feature"
    weight_normalization: str = "l1_abs"
    gene_scale: tuple[float, ...] | None = None


@dataclass(frozen=True)
class EncoderConfig:
    version: str; normalization: str; transformation: str; reference: str
    minimum_gene_coverage: float; minimum_cells: int; maximum_uncertainty: float
    required_modalities: tuple[str, ...]; optional_modalities: tuple[str, ...]; programs: tuple[ProgramSpec, ...]


class VersionedEncoder:
    def __init__(self, config: EncoderConfig):
        self.config = config
        payload = json.dumps(self.as_dict(include_hash=False), sort_keys=True, separators=(",", ":"))
        self.encoder_hash = hashlib.sha256(payload.encode()).hexdigest()

    @classmethod
    def from_yaml(cls, path: str | Path) -> "VersionedEncoder":
        raw = yaml.safe_load(Path(path).read_text())
        programs = tuple(ProgramSpec(p["name"], tuple(p["genes"]), tuple(p.get("weights", [1.0]*len(p["genes"]))), p.get("role", "candidate_feature"), p.get("weight_normalization", "l1_abs"), tuple(p["gene_scale"]) if "gene_scale" in p else None) for p in raw["programs"])
        return cls(EncoderConfig(programs=programs, required_modalities=tuple(raw.get("required_modalities", [])), optional_modalities=tuple(raw.get("optional_modalities", [])), **raw["encoder"]))

    def as_dict(self, include_hash=True):
        d = {"encoder": {k: getattr(self.config, k) for k in ("version","normalization","transformation","reference","minimum_gene_coverage","minimum_cells","maximum_uncertainty")},
             "required_modalities": list(self.config.required_modalities), "optional_modalities": list(self.config.optional_modalities),
             "programs": [{"name": p.name, "genes": list(p.genes), "weights": list(p.weights), "role": p.role, "weight_normalization": p.weight_normalization, "gene_scale": list(p.gene_scale) if p.gene_scale else None} for p in self.config.programs]}
        if include_hash: d["encoder_hash"] = self.encoder_hash
        return d

    def encode(self, expression: Mapping[str, float], *, sample_id: str, patient_id: str,
               cohort_id: str, timepoint: str, n_cells: int, available_modalities: Sequence[str],
               base_se: float = 1.0, input_scale: str = "log1p_CPM") -> LatentState:
        if input_scale != self.config.transformation:
            raise ValueError(f"encoder expects {self.config.transformation}, received {input_scale}; use encode_counts for raw counts")
        available = tuple(sorted(set(available_modalities)))
        effective_available = set(available)
        if effective_available.intersection({"scRNA-seq","snRNA-seq","bulk RNA-seq"}): effective_available.add("expression")
        missing = tuple(m for m in self.config.required_modalities if m not in effective_available)
        missing_optional = tuple(m for m in self.config.optional_modalities if m not in effective_available)
        modality_coverage = 1.0 if not self.config.required_modalities else sum(m in effective_available for m in self.config.required_modalities)/len(self.config.required_modalities)
        if set(available) - set(self.config.required_modalities) - set(self.config.optional_modalities):
            raise ValueError("incompatible modality supplied for this encoder")
        features = {}
        for p in self.config.programs:
            present = [i for i,g in enumerate(p.genes) if g in expression]
            coverage = len(present) / len(p.genes)
            if not present: value = float("nan")
            else:
                w = np.asarray([p.weights[i] for i in present], float)
                values = np.asarray([expression[p.genes[i]] for i in present], float)
                if p.weight_normalization == "l1_abs": w /= np.abs(w).sum()
                elif p.weight_normalization == "sum_present": w /= w.sum()
                elif p.weight_normalization != "none": raise ValueError(f"unknown weight normalization: {p.weight_normalization}")
                if p.gene_scale is not None: values /= np.asarray([p.gene_scale[i] for i in present])
                value = float(np.dot(w, values))
            parts = uncertainty_components(base_se,n_cells,coverage,modality_coverage,extrapolated=bool(missing))
            u = parts["uncertainty_score"]
            qc = "INSUFFICIENT" if n_cells < self.config.minimum_cells and coverage < self.config.minimum_gene_coverage else quality_label(u, self.config.maximum_uncertainty)
            if coverage < self.config.minimum_gene_coverage and qc == "PASS": qc = "LOW_CONFIDENCE"
            features[p.name] = LatentFeature(value, u, n_cells, coverage, modality_coverage, qc,
                parts["measurement_uncertainty"],parts["coverage_penalty"],parts["sampling_uncertainty"],parts["extrapolation_flag"])
        return LatentState(sample_id, patient_id, cohort_id, timepoint, features, available, missing_optional, missing, self.config.version, self.encoder_hash)

    def encode_counts(self, counts, **kwargs) -> LatentState:
        """Apply the configured raw-count transform before encoding."""
        if self.config.normalization != "library_size_CPM" or self.config.transformation != "log1p_CPM":
            raise ValueError("raw-count encoding requires normalization=library_size_CPM and transformation=log1p_CPM")
        raw = {str(k): float(v) for k,v in counts.items()}; total=sum(raw.values())
        if total <= 0: raise ValueError("raw counts have zero library size")
        transformed={k:float(np.log1p(v/total*1e6)) for k,v in raw.items()}
        return self.encode(transformed,input_scale="log1p_CPM",**kwargs)
