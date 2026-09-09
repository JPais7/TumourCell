"""Fail-closed sample-level temporal metadata validation."""
from dataclasses import dataclass
from collections import defaultdict

EVIDENCE=("DIRECT_SAMPLE_METADATA","DIRECT_STUDY_METADATA","SUPPLEMENTARY_TABLE","PUBLICATION_METHODS","PUBLICATION_TEXT_ONLY","INFERENCE","UNKNOWN")

@dataclass(frozen=True)
class LongitudinalSampleRecord:
    dataset_id: str
    patient_id: object = "UNKNOWN"
    sample_id: object = "UNKNOWN"
    biopsy_id: object = "UNKNOWN"
    timepoint_id: object = "UNKNOWN"
    timepoint_label: object = "UNKNOWN"
    biological_time_value: object = "UNKNOWN"
    biological_time_unit: object = "UNKNOWN"
    biological_time_anchor: object = "UNKNOWN"
    biological_time_source: object = "UNKNOWN"
    biological_time_evidence: str = "UNKNOWN"
    treatment: object = "UNKNOWN"
    treatment_start_available: object = "UNKNOWN"
    treatment_start_source: object = "UNKNOWN"
    sequencing_date: object = "UNKNOWN"
    run_date: object = "UNKNOWN"
    malignant_compartment: object = "UNKNOWN"
    clone_information: object = "UNKNOWN"
    record_source: str = "UNKNOWN"
    access_status: str = "UNKNOWN"

def _known(x): return x not in (None,"UNKNOWN","NOT_AVAILABLE","")

def validate_sample_identity(rows):
    mapping=defaultdict(set); samples=[]; biopsies=[]; duplicate_samples=[]; duplicate_biopsies=[]
    for r in rows:
        if _known(r.sample_id): mapping[r.sample_id].add(r.patient_id); samples.append(r.sample_id)
        if _known(r.biopsy_id): biopsies.append(r.biopsy_id)
    for s,p in mapping.items():
        if len(p)>1: duplicate_samples.append({'sample_id':s,'patients':sorted(map(str,p))})
    for s in set(samples):
        if samples.count(s)>1 and not any(x['sample_id']==s for x in duplicate_samples): duplicate_samples.append({'sample_id':s,'occurrences':samples.count(s)})
    for b in set(biopsies):
        if biopsies.count(b)>1: duplicate_biopsies.append({'biopsy_id':b,'occurrences':biopsies.count(b)})
    return {'sample_id_verified':bool(rows) and all(_known(r.sample_id) for r in rows),'patient_sample_mapping_verified':bool(rows) and all(_known(r.patient_id) for r in rows) and not any('patients' in x for x in duplicate_samples),'biopsy_identity_verified':bool(rows) and all(_known(r.biopsy_id) for r in rows),'duplicate_samples':duplicate_samples,'duplicate_biopsies':duplicate_biopsies}

def validate_temporal_values(rows):
    missing=[]; invalid=[]
    for r in rows:
        if not _known(r.biological_time_value) or not _known(r.biological_time_unit) or not _known(r.biological_time_anchor): missing.append(r.sample_id)
        elif not isinstance(r.biological_time_value,(int,float)) or isinstance(r.biological_time_value,bool): invalid.append(r.sample_id)
    return {'biological_time_observable':bool(rows) and not missing and not invalid,'missing_times':missing,'invalid_times':invalid}

def calculate_elapsed_intervals(rows):
    intervals=[]; missing=[]; invalid=[]; grouped=defaultdict(list)
    for r in rows: grouped[r.patient_id].append(r)
    for patient,rs in grouped.items():
        observed=[r.biological_time_value for r in rs]
        for a,b in zip(rs[:-1],rs[1:]):
            if isinstance(a.biological_time_value,(int,float)) and isinstance(b.biological_time_value,(int,float)) and b.biological_time_value < a.biological_time_value:
                invalid.append({'patient_id':patient,'from':a.sample_id,'to':b.sample_id,'delta_t':b.biological_time_value-a.biological_time_value})
        rs=sorted(rs,key=lambda r:r.biological_time_value if isinstance(r.biological_time_value,(int,float)) else float('inf'))
        for a,b in zip(rs[:-1],rs[1:]):
            if not isinstance(a.biological_time_value,(int,float)) or not isinstance(b.biological_time_value,(int,float)): missing.append({'patient_id':patient,'from':a.sample_id,'to':b.sample_id}); continue
            if a.biological_time_unit!=b.biological_time_unit or a.biological_time_anchor!=b.biological_time_anchor: invalid.append({'patient_id':patient,'from':a.sample_id,'to':b.sample_id,'reason':'unit_or_anchor_mismatch'}); continue
            d=b.biological_time_value-a.biological_time_value
            if d<=0: invalid.append({'patient_id':patient,'from':a.sample_id,'to':b.sample_id,'delta_t':d})
            else: intervals.append({'patient_id':patient,'from':a.sample_id,'to':b.sample_id,'delta_t':d,'unit':a.biological_time_unit})
    return {'intervals':intervals,'missing_intervals':missing,'invalid_intervals':invalid,'intervals_calculable':bool(intervals) and not missing,'intervals_valid':bool(intervals) and not invalid and not missing}

def validate_treatment_anchor(rows):
    return {'treatment_anchor_verified':bool(rows) and all(_known(r.treatment) and _known(r.biological_time_anchor) for r in rows)}

def classify_sample_level_temporal_observability(rows):
    identity=validate_sample_identity(rows); times=validate_temporal_values(rows); ints=calculate_elapsed_intervals(rows); treatment=validate_treatment_anchor(rows)
    repeated=any(sum(x.patient_id==r.patient_id for x in rows)>1 for r in rows if _known(r.patient_id)); malignant=bool(rows) and all(r.malignant_compartment is True for r in rows)
    if not rows: classification='UNKNOWN'
    elif repeated and identity['sample_id_verified'] and identity['patient_sample_mapping_verified'] and times['biological_time_observable'] and ints['intervals_valid'] and treatment['treatment_anchor_verified'] and malignant and not identity['duplicate_samples'] and not identity['duplicate_biopsies']: classification='CONTINUOUS_TIME_OBSERVABLE'
    elif repeated and all(_known(r.timepoint_label) for r in rows): classification='ORDERED_PHASE_ONLY'
    elif repeated: classification='LONGITUDINAL_BUT_INSUFFICIENT'
    else: classification='NOT_LONGITUDINAL'
    return {**identity,**times,**ints,**treatment,'repeated_patients':repeated,'malignant_compartment_verified':malignant,'classification':classification}

def classify_audit(rows):
    """Only validated sample collections are evidence; booleans cannot bypass validation."""
    return classify_sample_level_temporal_observability(rows) ['classification'] if isinstance(rows,(list,tuple)) else 'UNKNOWN'
