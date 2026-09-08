"""Unvalidated virtual perturbation hypotheses only."""
from dataclasses import dataclass
@dataclass(frozen=True)
class Perturbation:
    kind: str; target: str; direction: str; status: str="HYPOTHESIS"
    def __post_init__(self):
        if self.status != "HYPOTHESIS": raise ValueError("virtual perturbations must remain HYPOTHESIS")
