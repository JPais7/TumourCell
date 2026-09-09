"""Real-cohort validation using the frozen GSE246613 representation.

Temporal semantics are deliberately conservative: this cohort provides ordered
clinical phases, but no elapsed days.  Ordered phases are therefore suitable
for descriptive transitions and persistence, not for continuous-time dynamics.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ORDER={'Base':0,'PD1':1,'RTPD1':2}

def temporal_observations(npz_path):
    """Return an explicit, provenance-preserving temporal observation table.

    ``time_value`` is an ordinal clinical phase only.  It is never interpreted
    as a number of days; ``delta_t`` consequently remains ``None``.
    """
    x=np.load(npz_path,allow_pickle=True)
    out=[]
    for i,(patient,treatment) in enumerate(zip(x['patient'],x['treatment'])):
        patient=str(patient); treatment=str(treatment)
        phase=ORDER.get(treatment)
        out.append({'row_index':int(i),'patient_id':patient,
                    'sample_id':f'{patient}:observation_{i:03d}',
                    'timepoint_id':treatment,'time_value':phase,
                    'time_value_type':'ordered_clinical_phase',
                    'treatment':treatment,'state_row_index':int(i),
                    'response_label':str(x['response_group'][i]) if 'response_group' in x else None,
                    'temporal_status':'ORDERED_PHASE_ONLY'})
    return out

def temporal_transitions(observations):
    """Construct only adjacent, strictly forward patient transitions."""
    by_patient={}
    for row in observations: by_patient.setdefault(row['patient_id'],[]).append(row)
    transitions=[]; failures=[]
    for patient,rows in sorted(by_patient.items()):
        invalid=[r for r in rows if r.get('time_value') is None or r.get('timepoint_id') in (None,'')]
        if invalid:
            for r in invalid:
                failures.append({'patient_id':patient,'status':'MISSING_TIMEPOINT' if r.get('timepoint_id') in (None,'') else 'UNKNOWN_TIMEPOINT','timepoint_id':r.get('timepoint_id')})
            continue
        rows=sorted(rows,key=lambda r:(r['time_value'],r['row_index']))
        seen={}
        for r in rows:
            if r['time_value'] in seen:
                failures.append({'patient_id':patient,'status':'DUPLICATE_TIMEPOINT',
                                 'timepoint_id':r['timepoint_id']});
            seen[r['time_value']]=r
        if any(f['patient_id']==patient and f['status']=='DUPLICATE_TIMEPOINT' for f in failures):
            continue
        for history,target in zip(rows[:-1],rows[1:]):
            if target['time_value']<=history['time_value']:
                failures.append({'patient_id':patient,'status':'NON_FORWARD_OR_AMBIGUOUS'}); continue
            transitions.append({'patient_id':patient,
              'history_sample_id':history['sample_id'],'target_sample_id':target['sample_id'],
              'history_time':history['timepoint_id'],'target_time':target['timepoint_id'],
              'history_time_value':history['time_value'],'target_time_value':target['time_value'],
              'treatment_at_history':str(history['treatment']),'treatment_between':str(history['treatment'])+' -> '+str(target['treatment']),
              'delta_t':None,'delta_t_status':'NOT_AVAILABLE',
              'temporal_status':'ORDERED_PHASE_ONLY'})
    return transitions,failures

def _direction(value, epsilon=1e-8):
    return 'increase' if value>epsilon else ('decrease' if value<-epsilon else 'no_change')
def frozen_scores(npz_path, atlas_path, p8_path):
    x=np.load(npz_path,allow_pickle=True); genes=list(x['genes']); gi={g:i for i,g in enumerate(genes)}
    log=np.log1p(x['counts']/np.maximum(x['counts'].sum(1)[:,None],1)*1e6)
    atlas=json.loads(Path(atlas_path).read_text()); p8=json.loads(Path(p8_path).read_text())
    programs={s['state_id']:s['top_genes'] for s in atlas['states']}; programs['P8']=p8['genes']; names=list(programs)
    scores=np.zeros((len(log),len(names))); coverage={}
    for j,n in enumerate(names):
        present=[gi[g] for g in programs[n] if g in gi]; coverage[n]=len(present)/len(programs[n])
        if n=='P8':
            source=[i for i,g in enumerate(p8['genes']) if g in gi]; w=np.asarray(p8['weights'])[source]; w/=w.sum(); scale=np.asarray(p8['gene_scale'])[source]; scores[:,j]=(log[:,[gi[p8['genes'][i]] for i in source]]/scale)@w
        else: scores[:,j]=log[:,present].mean(1)
    return x,scores,names,coverage

def _metrics(pred,true,start,covs):
    err=pred-true; delta_true=true-start; delta_pred=pred-start; sd=np.asarray([np.sqrt(np.diag(c)) for c in covs])
    return {'MAE':float(np.mean(np.abs(err))),'RMSE':float(np.sqrt(np.mean(err**2))),'directional_accuracy':float(np.mean(np.sign(delta_true)==np.sign(delta_pred))),'directional_accuracy_per_feature':(np.sign(delta_true)==np.sign(delta_pred)).mean(0).tolist(),'predictive_interval_coverage_95':float(np.mean(np.abs(err)<=1.96*sd))}

def run_patient_held_out(npz_path,atlas_path,p8_path):
    x,z,names,coverage=frozen_scores(npz_path,atlas_path,p8_path)
    observations=temporal_observations(npz_path); transitions,failures=temporal_transitions(observations)
    rows=[]
    for tr in transitions:
        i=int(tr['history_sample_id'].split('_')[-1])
        j=int(tr['target_sample_id'].split('_')[-1])
        rows.append({**tr,'history':z[i],'observed':z[j],
                     'persistence_predicted':z[i],
                     'observed_delta_direction':[ _direction(v) for v in (z[j]-z[i]) ]})
    if not rows:
        return {'status':'NOT_ESTIMABLE','reason':'no eligible ordered transitions','coverage':coverage,'failures':failures}
    start=np.asarray([r['history'] for r in rows]); true=np.asarray([r['observed'] for r in rows])
    persistence=_metrics(start,true,start,[np.zeros((start.shape[1],start.shape[1])) for _ in rows])
    persistence['directional_accuracy']='NOT_APPLICABLE'
    persistence['directional_accuracy_per_feature']='NOT_APPLICABLE'
    persistence['predictive_interval_coverage_95']='NOT_ESTIMABLE'
    persistence['predicted_delta_direction']='no_change'
    return {'status':'ORDERED_PHASE_ONLY_NOT_CONTINUOUSLY_ESTIMABLE',
      'n_patients':len(set(r['patient_id'] for r in rows)),'n_observations':len(observations),
      'n_transitions':len(rows),'failures':failures,'coverage':coverage,
      'temporal_semantics':{'time_value_type':'ordered_clinical_phase','delta_t':'NOT_AVAILABLE',
                            'continuous_dynamics_fit':False,'duplicate_policy':'reject'},
      'persistence':persistence,
      'probabilistic_dynamics':{'status':'NOT_ESTIMABLE','reason':'exact elapsed time is unavailable'},
      'interval_coverage':{'50':'NOT_ESTIMABLE','80':'NOT_ESTIMABLE','95':'NOT_ESTIMABLE'},
      'response_used_as_predictor':False,
      'interpretation':'descriptive ordered-phase forecast only; no causal treatment interpretation'}

def forward_temporal_readiness(npz_path):
    """Assess, without fitting, whether a strict forward boundary has estimable treatments."""
    x=np.load(npz_path,allow_pickle=True); stage=np.asarray([ORDER.get(str(t),99) for t in x['treatment']]); transitions=[]
    for p in set(x['patient']):
        ids=[i for i,q in enumerate(x['patient']) if q==p]; ids.sort(key=lambda i:stage[i])
        transitions += [(stage[i],stage[j],str(x['treatment'][i])) for i,j in zip(ids[:-1],ids[1:])]
    # The first strict boundary has Base in history and PD1 as an unseen test treatment.
    return {'status':'NOT_ESTIMABLE','boundary_stage':1,'historical_transition_treatments':sorted(set(t for a,b,t in transitions if b<1)),'future_transition_treatments':sorted(set(t for a,b,t in transitions if a>=1)),'reason':'strict forward training at the first boundary cannot estimate the unseen PD1-conditioned transition without pooling or future leakage'}
