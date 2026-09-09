"""Explicit non-cascading digital-twin evidence gates."""
from enum import Enum
class GateStatus(str,Enum):
    PASSED="PASSED"; PARTIALLY_PASSED="PARTIALLY_PASSED"; NOT_PASSED="NOT_PASSED"; NOT_TESTED="NOT_TESTED"
GATES={"G1":"representation transport","G2":"external prediction","G3":"treatment-conditioned association","G4":"clone-resolved dynamics","G5":"spatial validation","G6":"counterfactual calibration","G7":"patient-informed state estimation","G8":"prospective prediction","G9":"synthetic model validation","G10":"observation-model validation"}
SIMULATION_PREREQUISITE_GATES=("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10")
def evaluate_gates(evidence):
    result={g:GateStatus(evidence.get(g,GateStatus.NOT_TESTED)) for g in GATES}
    simulation_enabled=simulation_allowed({g:result[g].value for g in GATES})
    return {"gates":{g:{"name":GATES[g],"status":result[g].value} for g in GATES},"simulation_enabled":simulation_enabled,
            "clinical_recommendations_enabled":False}

def simulation_allowed(evidence=None):
    """Always fail closed until every independent gate is explicitly passed."""
    evidence=evidence or {}
    return all(evidence.get(g)==GateStatus.PASSED.value for g in SIMULATION_PREREQUISITE_GATES) and False
