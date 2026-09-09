"""Known-truth generators for selection, plasticity and confounding tests."""
import numpy as np
def generate(mechanism,n_patients=40,seed=0):
    allowed={"selection","plasticity","selection_plus_plasticity","sampling_depth","spatial_environment"}
    if mechanism not in allowed: raise ValueError(f"unknown mechanism: {mechanism}")
    rng=np.random.default_rng(seed); clone=rng.binomial(1,.5,n_patients); state=rng.normal(size=n_patients); env=rng.normal(size=n_patients)
    follow=state.copy(); clone_follow=clone.copy()
    if mechanism in {"selection","selection_plus_plasticity"}: clone_follow=np.where(rng.random(n_patients)<.7,1,clone)
    if mechanism in {"plasticity","selection_plus_plasticity"}: follow+=1
    if mechanism=="sampling_depth": follow+=rng.normal(scale=1/np.sqrt(rng.integers(5,500,n_patients)))
    if mechanism=="spatial_environment": follow+=env
    return {"baseline_state":state,"followup_state":follow,"baseline_clone":clone,"followup_clone":clone_follow,"environment":env,"truth":mechanism}

def clonal_selection(seed=0):
    rng=np.random.default_rng(seed); baseline=np.array([0]*70+[1]*30); follow=np.array([0]*20+[1]*80)
    return {"baseline_clone":rng.permutation(baseline),"followup_clone":rng.permutation(follow),"truth":"selection"}
def transcriptional_plasticity(seed=0):
    rng=np.random.default_rng(seed); clones=rng.integers(0,2,100)
    return {"baseline_clone":clones,"followup_clone":clones.copy(),"baseline_state":rng.normal(size=100),"followup_state":rng.normal(size=100)+1,"truth":"plasticity"}
def selection_plus_plasticity(seed=0):
    d=clonal_selection(seed); rng=np.random.default_rng(seed); d.update(baseline_state=rng.normal(size=100),followup_state=rng.normal(size=100)+1); d["truth"]="selection_plus_plasticity"; return d
def recover(data):
    clone_change=float(np.mean(data["baseline_clone"]!=data["followup_clone"]))
    state_change=float(np.mean(data["followup_state"]-data["baseline_state"]))
    delta=data["followup_state"]-data["baseline_state"]
    env_corr=float(np.corrcoef(data["environment"],delta)[0,1]) if np.std(delta)>0 else 0.0
    if abs(env_corr)>.8:return "spatial_environment"
    if clone_change>.1 and state_change>.5:return "selection_plus_plasticity"
    if clone_change>.1:return "selection"
    if state_change>.5:return "plasticity"
    return "sampling_depth"
