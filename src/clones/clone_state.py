"""Separate selection, plasticity and environmental effects."""
from dataclasses import dataclass
NOT_IDENTIFIABLE="NOT_IDENTIFIABLE"
@dataclass(frozen=True)
class CloneState:
    clone: str | None; cellular_state: tuple[float,...]; environment: tuple[float,...]
def identify_mechanisms(*, has_clone_data, assignment_quality=0.0, clones_by_timepoint=None,
                        minimum_quality=0.8, minimum_abundance=0.05, temporal_overlap=False,
                        sampling_qc=False, environment_measured=False):
    """Return estimability only when clone identity is observed and comparable."""
    clones_by_timepoint = clones_by_timepoint or {}
    timepoints=list(clones_by_timepoint)
    abundant = all(any(v >= minimum_abundance for v in d.values()) for d in clones_by_timepoint.values()) if timepoints else False
    comparable = len(timepoints) >= 2 and len(set.intersection(*(set(d) for d in clones_by_timepoint.values()))) > 0 if len(timepoints)>=2 else False
    enough = has_clone_data and assignment_quality >= minimum_quality and abundant and comparable and temporal_overlap and sampling_qc
    if not enough:
        return {"selection_status":NOT_IDENTIFIABLE,"plasticity_status":NOT_IDENTIFIABLE,
                "environmental_effect_status":"ESTIMABLE" if environment_measured else NOT_IDENTIFIABLE}
    return {"selection_status":"ESTIMABLE","plasticity_status":"ESTIMABLE" if comparable and temporal_overlap else NOT_IDENTIFIABLE,
            "environmental_effect_status":"ESTIMABLE" if environment_measured else NOT_IDENTIFIABLE}

def classify_observation(*, expression_before_after=False, clone_evidence=False, environment_evidence=False):
    if clone_evidence: return "selection"
    if environment_evidence: return "environmental_shift"
    if expression_before_after: return "observed_state_shift"
    return "NOT_IDENTIFIABLE"
