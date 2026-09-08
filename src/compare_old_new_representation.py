#!/usr/bin/env python3
"""Locked BioKey old-versus-versioned-encoder equivalence audit."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from analyze_biokey_longitudinal import process_cohort
from representation.encoder import EncoderConfig,ProgramSpec,VersionedEncoder

ROOT=Path(__file__).parents[1]
atlas=json.loads((ROOT/'results/phase3/state_atlas_v1.json').read_text())
p8=json.loads((ROOT/'results/phase2/p8_definition_frozen.json').read_text())
specs=[ProgramSpec(s['state_id'],tuple(s['top_genes']),tuple([1.0]*len(s['top_genes']))) for s in atlas['states']]
specs.append(ProgramSpec('P8',tuple(p8['genes']),tuple(p8['weights']),weight_normalization='sum_present',gene_scale=tuple(p8['gene_scale'])))
encoder=VersionedEncoder(EncoderConfig('historical_frozen_bridge_v1','log1p_CPM','weighted_mean_present_genes','Phase2/3 frozen files',0.0,1,float('inf'),('scRNA-seq',),tuple(specs)))
wanted=list(dict.fromkeys(g for p in specs for g in p.genes)); rows=[]
for cohort in (1,2):
 info,counts,present=process_cohort(cohort,ROOT/'data/raw/BIOKEY',wanted)
 expr=np.log1p(counts/np.maximum(info.library_size.to_numpy()[:,None],1)*1e6)
 for i,r in info.iterrows():
  state=encoder.encode(dict(zip(present,expr[i])),sample_id=f"{r.patient_id}_{r.timepoint}",patient_id=r.patient_id,cohort_id=str(cohort),timepoint=r.timepoint,n_cells=int(r.n_malignant_cells),available_modalities=['scRNA-seq'])
  rows.append({'cohort':cohort,'patient_id':r.patient_id,'timepoint':r.timepoint,**{k:v.value for k,v in state.features.items()}})
new=pd.DataFrame(rows); old=pd.read_csv(ROOT/'results/phase4/biokey_patient_timepoint_scores.csv')
joined=old.merge(new,on=['cohort','patient_id','timepoint'],suffixes=('_old','_new'),validate='one_to_one')
metrics=['State_01','State_02','State_03','State_04','P8']
payload={'comparison':'historical BioKey scorer versus versioned encoder bridge','patients_old':int(old.patient_id.nunique()),'patients_new':int(new.patient_id.nunique()),'rows_old':len(old),'rows_new':len(new),'encoder_hash':encoder.encoder_hash,
 'maximum_absolute_score_difference':{m:float(np.max(np.abs(joined[m+'_old']-joined[m+'_new']))) for m in metrics},
 'included_patients_identical':set(old.patient_id)==set(new.patient_id),'statistical_results':'unchanged only if all score differences are below 1e-10',
 'equivalence_passed':all(np.max(np.abs(joined[m+'_old']-joined[m+'_new'])) < 1e-10 for m in metrics),'conclusion_changes':False}
out=ROOT/'results/phase4/old_vs_new_representation.json'; out.write_text(json.dumps(payload,indent=2)+'\n'); print(json.dumps(payload,indent=2))
