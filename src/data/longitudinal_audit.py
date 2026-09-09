"""Sample-level metadata audit primitives for Phase 4.5.1."""
from dataclasses import dataclass, asdict

EVIDENCE = ("DIRECT_SAMPLE_METADATA", "DIRECT_STUDY_METADATA", "SUPPLEMENTARY_TABLE",
            "PUBLICATION_METHODS", "PUBLICATION_TEXT_ONLY", "INFERENCE", "UNKNOWN")

@dataclass(frozen=True)
class LongitudinalMetadataAudit:
    dataset_id: str
    patient_id: object = "UNKNOWN"
    sample_id: object = "UNKNOWN"
    repeated_samples: object = "UNKNOWN"
    elapsed_time_observable: object = "UNKNOWN"
    elapsed_time_source: str = "UNKNOWN"
    elapsed_time_unit: str = "UNKNOWN"
    elapsed_time_anchor: str = "UNKNOWN"
    timepoint_order_observable: object = "UNKNOWN"
    biological_time_distinguishable_from_sequencing_time: object = "UNKNOWN"
    treatment: object = "UNKNOWN"
    malignant_compartment: object = "UNKNOWN"
    clone_information: object = "UNKNOWN"
    evidence: tuple = ()
    classification: str = "UNKNOWN"
    confidence: str = "LOW"
    limitations: tuple = ()
    next_action: str = ""

def classify_audit(a: LongitudinalMetadataAudit):
    """Fail-closed classification; publication phase text cannot create elapsed time."""
    if a.repeated_samples is False: return "NOT_LONGITUDINAL"
    if a.repeated_samples is not True: return "UNKNOWN"
    if (a.elapsed_time_observable is True and a.patient_id is True and
        a.sample_id is True and a.treatment is True and a.malignant_compartment is True):
        return "CONTINUOUS_TIME_OBSERVABLE"
    if a.timepoint_order_observable is True: return "ORDERED_PHASE_ONLY"
    return "LONGITUDINAL_BUT_INSUFFICIENT"

def audit_record(a):
    d=asdict(a); d['evidence']=list(a.evidence); d['limitations']=list(a.limitations); return d
