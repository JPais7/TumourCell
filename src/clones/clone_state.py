"""Separate selection, plasticity and environmental effects."""
from dataclasses import dataclass
NOT_IDENTIFIABLE="NOT_IDENTIFIABLE"
from enum import Enum
class MechanismStatus(str,Enum): SUPPORTED="SUPPORTED"; NOT_SUPPORTED="NOT_SUPPORTED"; INDETERMINATE="INDETERMINATE"; NOT_IDENTIFIABLE="NOT_IDENTIFIABLE"
@dataclass(frozen=True)
class ObservedShiftAssessment:
    observed_state_shift: str; composition_shift: str; selection: str; plasticity: str
    environmental_shift: str; sampling_confounding: str; technical_confounding: str
    identifiability_status: str; evidence: dict; reasons: tuple[str,...]
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
    """Deprecated fail-closed wrapper; booleans never infer mechanisms."""
    import warnings
    warnings.warn("classify_observation is deprecated; use assess_observed_shift", DeprecationWarning, stacklevel=2)
    return assess_observed_shift(expression_before_after=expression_before_after, clone_evidence=clone_evidence, environment_evidence=environment_evidence)

def assess_observed_shift(*, expression_before_after=False, clone_evidence=False, assignment_quality=0.0,
                          temporal_overlap=False, abundance_sufficient=False, sampling_comparable=False,
                          clone_frequency_changed=False, within_clone_state_change=False,
                          environment_measured=False, sampling_problem=False, technical_problem=False):
    reasons=[]; evidence={"expression_before_after":expression_before_after,"clone_evidence":clone_evidence}
    selection = MechanismStatus.SUPPORTED.value if (clone_evidence and assignment_quality>=.8 and temporal_overlap and abundance_sufficient and sampling_comparable and clone_frequency_changed) else (MechanismStatus.INDETERMINATE.value if clone_evidence else MechanismStatus.NOT_IDENTIFIABLE.value)
    plasticity = MechanismStatus.SUPPORTED.value if (clone_evidence and temporal_overlap and within_clone_state_change and sampling_comparable) else (MechanismStatus.INDETERMINATE.value if clone_evidence else MechanismStatus.NOT_IDENTIFIABLE.value)
    environmental = MechanismStatus.SUPPORTED.value if environment_measured else MechanismStatus.NOT_IDENTIFIABLE.value
    sampling = MechanismStatus.SUPPORTED.value if sampling_problem else MechanismStatus.NOT_SUPPORTED.value
    technical = MechanismStatus.SUPPORTED.value if technical_problem else MechanismStatus.NOT_SUPPORTED.value
    state = MechanismStatus.SUPPORTED.value if expression_before_after else MechanismStatus.NOT_IDENTIFIABLE.value
    composition = MechanismStatus.INDETERMINATE.value if expression_before_after and not clone_evidence else MechanismStatus.NOT_IDENTIFIABLE.value
    identifiable=[x for x in (selection,plasticity,environmental,sampling,technical) if x==MechanismStatus.SUPPORTED.value]
    status=MechanismStatus.INDETERMINATE.value if len(identifiable)>1 else (MechanismStatus.SUPPORTED.value if identifiable else MechanismStatus.NOT_IDENTIFIABLE.value)
    if not clone_evidence: reasons.append("no clone identity: selection/plasticity cannot be separated")
    if expression_before_after: reasons.append("expression change supports observed_state_shift only")
    return ObservedShiftAssessment(state,composition,selection,plasticity,environmental,sampling,technical,status,evidence,tuple(reasons))
