"""Conservative qualification of longitudinal cancer datasets.

This module classifies observability; it does not fit a temporal model.
Unknown metadata stays unknown and sequencing dates are never biological time.
"""
from dataclasses import dataclass, asdict
from typing import Any, Iterable

CLASSES = ("CONTINUOUS_TIME_OBSERVABLE", "ORDERED_PHASE_ONLY",
           "LONGITUDINAL_BUT_INSUFFICIENT", "NOT_LONGITUDINAL", "UNKNOWN")

@dataclass(frozen=True)
class TemporalDatasetQualification:
    dataset_id: str
    disease: str = "UNKNOWN"
    subtype: str = "UNKNOWN"
    assay: str = "UNKNOWN"
    single_cell: Any = "UNKNOWN"
    single_nucleus: Any = "UNKNOWN"
    bulk: Any = "UNKNOWN"
    spatial: Any = "UNKNOWN"
    patient_id_available: Any = "UNKNOWN"
    sample_id_available: Any = "UNKNOWN"
    repeated_patient_samples: Any = "UNKNOWN"
    timestamp_available: Any = "UNKNOWN"
    date_available: Any = "UNKNOWN"
    days_from_treatment_available: Any = "UNKNOWN"
    elapsed_time_observable: Any = "UNKNOWN"
    ordered_phase_available: Any = "UNKNOWN"
    treatment_metadata_available: Any = "UNKNOWN"
    malignant_compartment_available: Any = "UNKNOWN"
    malignant_cells_available: Any = "UNKNOWN"
    clone_information_available: Any = "UNKNOWN"
    spatial_information_available: Any = "UNKNOWN"
    qualification: str = "UNKNOWN"
    evidence: str = ""
    limitations: str = ""
    provenance: str = ""

    def __post_init__(self):
        if self.qualification not in CLASSES:
            raise ValueError(f"unknown qualification: {self.qualification}")

def qualify_temporal_dataset(*, patient_id_available, sample_id_available,
    repeated_patient_samples, elapsed_time_observable, ordered_phase_available,
    treatment_metadata_available, malignant_compartment_available,
    data_available=True, **kwargs):
    """Deterministic, fail-closed qualification rule."""
    if not data_available:
        cls = "UNKNOWN"
    elif repeated_patient_samples is False:
        cls = "NOT_LONGITUDINAL"
    elif repeated_patient_samples is True and elapsed_time_observable is True \
         and patient_id_available is True and sample_id_available is True \
         and treatment_metadata_available is True \
         and malignant_compartment_available is True:
        cls = "CONTINUOUS_TIME_OBSERVABLE"
    elif repeated_patient_samples is True and ordered_phase_available is True:
        cls = "ORDERED_PHASE_ONLY"
    elif repeated_patient_samples is True:
        cls = "LONGITUDINAL_BUT_INSUFFICIENT"
    else:
        cls = "UNKNOWN"
    return cls

def temporal_dataset_gate(q: TemporalDatasetQualification):
    """Return explicit gate fields; no longitudinal shortcut can pass."""
    return {"continuous_time_available": q.qualification == "CONTINUOUS_TIME_OBSERVABLE",
            "ordered_phase_available": q.ordered_phase_available,
            "repeated_patient_samples": q.repeated_patient_samples,
            "treatment_observable": q.treatment_metadata_available,
            "malignant_compartment_observable": q.malignant_compartment_available,
            "sufficient_for_forecasting": q.qualification == "CONTINUOUS_TIME_OBSERVABLE",
            "classification": q.qualification}

def qualification_record(q: TemporalDatasetQualification):
    d = asdict(q); d["classification"] = d.pop("qualification"); return d

def rank_qualification(q: TemporalDatasetQualification):
    """Transparent qualitative ranking, not an arbitrary numeric score."""
    if q.qualification == "CONTINUOUS_TIME_OBSERVABLE" and q.malignant_compartment_available is True:
        return "BEST_CANDIDATE"
    if q.qualification == "ORDERED_PHASE_ONLY": return "PROMISING"
    if q.qualification == "LONGITUDINAL_BUT_INSUFFICIENT": return "PARTIAL"
    return "INSUFFICIENT"
