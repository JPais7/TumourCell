"""Group/temporal validation for probabilistic transitions."""
from __future__ import annotations
import numpy as np
def validate_transition_model(model, current_state, next_state, treatment, delta_t, patient_id):
    z=np.asarray(current_state,float); y=np.asarray(next_state,float); pid=np.asarray(patient_id)
    predictions=[]; truths=[]; starts=[]; rows=[]; ll=[]; covs=[]
    for patient in dict.fromkeys(pid):
        train=pid!=patient
        if train.sum()<5: continue
        fitted=type(model)().fit(z[train],y[train],np.asarray(treatment)[train],np.asarray(delta_t)[train],pid[train])
        idx=np.flatnonzero(~train)
        for i in idx:
            dist=fitted.predict_distribution(z[i],treatment[i],float(delta_t[i])); pred=dist['mean'][0]
            predictions.append(pred); truths.append(y[i]); starts.append(z[i]); covs.append(dist['covariance']); ll.append(fitted.log_likelihood(z[i],y[i],treatment[i],float(delta_t[i])))
            rows.append({'patient_id':str(patient),'treatment':str(treatment[i]),'delta_t':float(delta_t[i]),'predicted_mean':pred.tolist(),'observed_next_state':y[i].tolist(),'predictive_covariance':dist['covariance'].tolist(),'residual':(y[i]-pred).tolist(),'fold_id':str(patient)})
    if not predictions: raise ValueError("not enough patients for leave-one-patient-out validation")
    p=np.asarray(predictions); truth=np.asarray(truths); starts=np.asarray(starts); err=p-truth; rmse=float(np.sqrt(np.mean(err**2))); mae=float(np.mean(np.abs(err)))
    sd=np.asarray([np.sqrt(np.diag(c)) for c in covs]); coverage=float(np.mean(np.abs(err)<=1.96*sd))
    true_delta=truth-starts; pred_delta=p-starts; directional=np.sign(pred_delta)==np.sign(true_delta)
    per_feature=directional.mean(axis=0).tolist()
    return {"n_held_out_patients":len(set(r['patient_id'] for r in rows)),"n_held_out_transitions":len(p),"fold_count":len(rows),"MAE":mae,"RMSE":rmse,"mean_predictive_log_likelihood":float(np.mean(ll)),"predictive_interval_coverage_95":coverage,"directional_accuracy":float(directional.mean()),"directional_accuracy_per_feature":per_feature,"zero_change_rule":"zero/zero counts as correct; zero/nonzero as incorrect","split":"leave_one_patient_out","unit":"patient","fold_predictions":rows}
