"""Group/temporal validation for probabilistic transitions."""
from __future__ import annotations
import numpy as np
def validate_transition_model(model, current_state, next_state, treatment, delta_t, patient_id):
    z=np.asarray(current_state,float); y=np.asarray(next_state,float); pid=np.asarray(patient_id)
    predictions=[]; truths=[]; ll=[]
    for patient in dict.fromkeys(pid):
        train=pid!=patient
        if train.sum()<5: continue
        fitted=type(model)().fit(z[train],y[train],np.asarray(treatment)[train],np.asarray(delta_t)[train],pid[train])
        idx=np.flatnonzero(~train)
        for i in idx:
            dist=fitted.predict_distribution(z[i],treatment[i],float(delta_t[i])); predictions.append(dist['mean'][0]); truths.append(y[i]); ll.append(fitted.log_likelihood(z[i],y[i],treatment[i],float(delta_t[i])))
    if not predictions: raise ValueError("not enough patients for leave-one-patient-out validation")
    p=np.asarray(predictions); truth=np.asarray(truths); err=p-truth; rmse=float(np.sqrt(np.mean(err**2))); mae=float(np.mean(np.abs(err)))
    sd=np.asarray([np.sqrt(np.diag(model.residual_covariance[str(t)])) for t in np.asarray(treatment)[-len(p):]]) if model.residual_covariance else np.ones_like(p)
    return {"n_held_out_transitions":len(p),"MAE":mae,"RMSE":rmse,"predictive_log_likelihood":float(np.mean(ll)),"prediction_interval_coverage_95":float(np.mean(np.abs(err)<=1.96*sd)),"directional_accuracy":float(np.mean(np.sign(p)==np.sign(truth))),"split":"leave_one_patient_out","unit":"patient"}
