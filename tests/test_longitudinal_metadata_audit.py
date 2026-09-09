from src.data.longitudinal_audit import LongitudinalSampleRecord as R, classify_sample_level_temporal_observability as classify, calculate_elapsed_intervals

def row(p='P1',s='S1',b='B1',t=0,label='baseline',**kw):
    d=dict(dataset_id='x',patient_id=p,sample_id=s,biopsy_id=b,timepoint_id=s,timepoint_label=label,biological_time_value=t,biological_time_unit='days',biological_time_anchor='treatment_start',treatment='NAC',malignant_compartment=True)
    d.update(kw); return R(**d)

def test_real_days_allow_continuous_time(): assert classify([row(),row(s='S2',b='B2',t=7,label='on')])['classification']=='CONTINUOUS_TIME_OBSERVABLE'
def test_duplicate_sample_rows_are_detected(): assert classify([row(),row(s='S1',b='B2',t=10)])['duplicate_samples']
def test_sample_conflicting_patients_fails_mapping(): assert classify([row(),row(p='P2',s='S1',b='B2',t=10)])['patient_sample_mapping_verified'] is False
def test_duplicate_biopsy_is_detected(): assert classify([row(),row(s='S2',b='B1',t=10)])['duplicate_biopsies']
def test_same_time_samples_are_not_adjacent_valid_intervals(): assert classify([row(),row(s='S2',b='B2',t=0,label='rep')])['intervals_valid'] is False
def test_legitimate_irregular_intervals_preserved(): assert [x['delta_t'] for x in calculate_elapsed_intervals([row(),row(s='S2',b='B2',t=7),row(s='S3',b='B3',t=21)])['intervals']]==[7,14]
def test_negative_interval_rejected(): assert classify([row(t=21),row(s='S2',b='B2',t=7)])['invalid_intervals']
def test_missing_time_not_interpolated(): assert classify([row(),row(s='S2',b='B2',t='UNKNOWN'),row(s='S3',b='B3',t=21)])['missing_times']
def test_sequencing_date_does_not_count(): assert classify([row(sequencing_date='2024-01-01'),row(s='S2',b='B2',t='UNKNOWN',sequencing_date='2024-01-20')])['classification']!='CONTINUOUS_TIME_OBSERVABLE'
def test_missing_anchor_blocks_pass(): assert classify([row(biological_time_anchor='UNKNOWN'),row(s='S2',b='B2',t=7)])['classification']!='CONTINUOUS_TIME_OBSERVABLE'
def test_cycle_label_does_not_create_time(): assert classify([row(t='UNKNOWN',biological_time_unit='UNKNOWN',biological_time_anchor='UNKNOWN',timepoint_label='cycle1'),row(s='S2',b='B2',t='UNKNOWN',biological_time_unit='UNKNOWN',biological_time_anchor='UNKNOWN',timepoint_label='cycle2')])['classification']=='ORDERED_PHASE_ONLY'
def test_unknown_records_fail_closed(): assert classify([])['classification']=='UNKNOWN'
def test_phase_only_is_not_continuous(): assert classify([row(t='UNKNOWN',biological_time_unit='UNKNOWN',biological_time_anchor='UNKNOWN'),row(s='S2',b='B2',t='UNKNOWN',biological_time_unit='UNKNOWN',biological_time_anchor='UNKNOWN',label='on')])['classification']=='ORDERED_PHASE_ONLY'
def test_future_row_does_not_change_earlier_value():
    early=[row(),row(s='S2',b='B2',t=7)]; assert classify(early)['intervals'][0]['delta_t']==7
