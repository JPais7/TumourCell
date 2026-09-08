"""Separate selection, plasticity and environmental effects."""
from dataclasses import dataclass
NOT_IDENTIFIABLE="NOT_IDENTIFIABLE"
@dataclass(frozen=True)
class CloneState:
    clone: str | None; cellular_state: tuple[float,...]; environment: tuple[float,...]
def identify_mechanisms(*, has_clone_data, repeated_clones=False, environment_measured=False):
    if not has_clone_data:
        return {"selection_status":NOT_IDENTIFIABLE,"plasticity_status":NOT_IDENTIFIABLE,
                "environmental_effect_status":"ESTIMABLE" if environment_measured else NOT_IDENTIFIABLE}
    return {"selection_status":"ESTIMABLE","plasticity_status":"ESTIMABLE" if repeated_clones else NOT_IDENTIFIABLE,
            "environmental_effect_status":"ESTIMABLE" if environment_measured else NOT_IDENTIFIABLE}
