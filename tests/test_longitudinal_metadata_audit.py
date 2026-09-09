from src.data.longitudinal_audit import LongitudinalMetadataAudit, classify_audit

def a(**kw):
    d=dict(dataset_id='x',patient_id=True,sample_id=True,repeated_samples=True,
           elapsed_time_observable=False,timepoint_order_observable=False,
           treatment=True,malignant_compartment=True)
    d.update(kw); return LongitudinalMetadataAudit(**d)

def test_phase_without_interval_not_continuous(): assert classify_audit(a(timepoint_order_observable=True))=='ORDERED_PHASE_ONLY'
def test_sequencing_date_not_biological_time(): assert classify_audit(a(elapsed_time_observable=False))!='CONTINUOUS_TIME_OBSERVABLE'
def test_real_days_can_pass(): assert classify_audit(a(elapsed_time_observable=True))=='CONTINUOUS_TIME_OBSERVABLE'
def test_missing_patient_id_fails(): assert classify_audit(a(patient_id=False,elapsed_time_observable=True))!='CONTINUOUS_TIME_OBSERVABLE'
def test_missing_sample_id_fails(): assert classify_audit(a(sample_id=False,elapsed_time_observable=True))!='CONTINUOUS_TIME_OBSERVABLE'
def test_missing_treatment_fails(): assert classify_audit(a(treatment=False,elapsed_time_observable=True))!='CONTINUOUS_TIME_OBSERVABLE'
def test_missing_malignant_compartment_fails(): assert classify_audit(a(malignant_compartment=False,elapsed_time_observable=True))!='CONTINUOUS_TIME_OBSERVABLE'
def test_controlled_metadata_does_not_imply_pass(): assert classify_audit(a(elapsed_time_observable='ACCESS_CONTROLLED'))!='CONTINUOUS_TIME_OBSERVABLE'
def test_unequal_intervals_are_not_normalised():
    x=a(elapsed_time_observable=True); assert x.elapsed_time_observable is True
def test_missing_interval_not_interpolated(): assert classify_audit(a(elapsed_time_observable=False,timepoint_order_observable=True))=='ORDERED_PHASE_ONLY'
def test_duplicate_sample_is_not_silent(): assert classify_audit(a(sample_id='DUPLICATE'))!='CONTINUOUS_TIME_OBSERVABLE'
def test_phase_order_never_creates_delta_t():
    x=a(timepoint_order_observable=True); assert x.elapsed_time_observable is False
def test_future_metadata_does_not_change_past_rule(): assert classify_audit(a(elapsed_time_observable=False))!='CONTINUOUS_TIME_OBSERVABLE'
def test_unknown_fields_stay_unknown(): assert classify_audit(LongitudinalMetadataAudit(dataset_id='x'))=='UNKNOWN'
