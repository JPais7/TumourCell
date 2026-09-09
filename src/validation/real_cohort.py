"""Minimal real-cohort benchmark using frozen GSE246613 pseudobulks."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from src.dynamics.discrete import ProbabilisticTransition

ORDER={'Base':0,'PD1':1,'RTPD1':2}
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
    x,z,names,coverage=frozen_scores(npz_path,atlas_path,p8_path); rows=[]; failures=[]
    for patient in sorted(set(x['patient'])):
        idx=[i for i,p in enumerate(x['patient']) if p==patient]; idx.sort(key=lambda i:ORDER.get(str(x['treatment'][i]),99))
        train_idx=[i for i,p in enumerate(x['patient']) if p!=patient]; transitions=[]
        for i,j in zip(idx[:-1],idx[1:]):
            if ORDER.get(str(x['treatment'][i]),99)>=ORDER.get(str(x['treatment'][j]),99): continue
            transitions.append((i,j))
        if not transitions: continue
        tr=[]
        for i in train_idx:
            for j in train_idx:
                if x['patient'][i]!=x['patient'][j] or ORDER.get(str(x['treatment'][i]),99)+1!=ORDER.get(str(x['treatment'][j]),99): continue
                tr.append((i,j))
        if not tr: failures.append({'patient_id':str(patient),'reason':'no_training_transitions'}); continue
        z0=np.asarray([z[i] for i,j in tr]); z1=np.asarray([z[j] for i,j in tr]); tt=np.asarray([str(x['treatment'][i]) for i,j in tr]); dt=np.ones(len(tr)); pid=np.asarray([str(x['patient'][i]) for i,j in tr])
        if any(sum(tt==t)<3 for t in set(tt)): failures.append({'patient_id':str(patient),'reason':'treatment_transition_count_below_3'}); continue
        try: model=ProbabilisticTransition().fit(z0,z1,tt,dt,pid)
        except ValueError as e: failures.append({'patient_id':str(patient),'reason':str(e)}); continue
        for i,j in transitions:
            treatment=str(x['treatment'][i]);
            if treatment not in model.treatment_levels: failures.append({'patient_id':str(patient),'reason':'unseen_treatment'}); continue
            dist=model.predict_distribution(z[i],treatment,1.0); pred=dist['mean'][0]
            rows.append({'patient_id':str(patient),'treatment':treatment,'history_time':str(x['treatment'][i]),'target_time':str(x['treatment'][j]),'predicted':pred,'observed':z[j],'history':z[i],'covariance':dist['covariance']})
    if not rows: return {'status':'NOT_ESTIMABLE','reason':'no eligible patient-held-out transitions','coverage':coverage,'failures':failures}
    pred=np.asarray([r['predicted'] for r in rows]); true=np.asarray([r['observed'] for r in rows]); start=np.asarray([r['history'] for r in rows]); cov=[r['covariance'] for r in rows]
    persistence=_metrics(start,true,start,cov); dynamics=_metrics(pred,true,start,cov)
    return {'status':'SYNTHETICALLY_UNVALIDATED_REAL_BENCHMARK','n_patients':len(set(r['patient_id'] for r in rows)),'n_transitions':len(rows),'failures':failures,'coverage':coverage,'persistence':persistence,'probabilistic_dynamics':dynamics,'treatment_counts':{t:sum(r['treatment']==t for r in rows) for t in set(r['treatment'] for r in rows)},'response_used_as_predictor':False,'interpretation':'predictive benchmark only; no causal treatment interpretation'}

def forward_temporal_readiness(npz_path):
    """Assess, without fitting, whether a strict forward boundary has estimable treatments."""
    x=np.load(npz_path,allow_pickle=True); stage=np.asarray([ORDER.get(str(t),99) for t in x['treatment']]); transitions=[]
    for p in set(x['patient']):
        ids=[i for i,q in enumerate(x['patient']) if q==p]; ids.sort(key=lambda i:stage[i])
        transitions += [(stage[i],stage[j],str(x['treatment'][i])) for i,j in zip(ids[:-1],ids[1:])]
    # The first strict boundary has Base in history and PD1 as an unseen test treatment.
    return {'status':'NOT_ESTIMABLE','boundary_stage':1,'historical_transition_treatments':sorted(set(t for a,b,t in transitions if b<1)),'future_transition_treatments':sorted(set(t for a,b,t in transitions if a>=1)),'reason':'strict forward training at the first boundary cannot estimate the unseen PD1-conditioned transition without pooling or future leakage'}
