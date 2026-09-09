from src.data.temporal_qualification import *

def make(**kw):
    base=dict(dataset_id='x', patient_id_available=True, sample_id_available=True,
              repeated_patient_samples=True, timestamp_available=False,
              date_available=False, days_from_treatment_available=False,
              elapsed_time_observable=False, ordered_phase_available=False,
              treatment_metadata_available=True, malignant_compartment_available=True)
    base.update(kw); rule=dict(base); data_available=rule.pop('data_available', True)
    return TemporalDatasetQualification(**rule, qualification=qualify_temporal_dataset(**rule, data_available=data_available))

def test_real_elapsed_time_passes():
    assert make(elapsed_time_observable=True).qualification == 'CONTINUOUS_TIME_OBSERVABLE'

def test_phases_only_are_not_continuous():
    q=make(ordered_phase_available=True); assert q.qualification == 'ORDERED_PHASE_ONLY'

def test_one_sample_is_not_longitudinal():
    assert make(repeated_patient_samples=False).qualification == 'NOT_LONGITUDINAL'

def test_repeated_without_patient_id_is_insufficient():
    assert make(patient_id_available=False).qualification == 'LONGITUDINAL_BUT_INSUFFICIENT'

def test_sequencing_date_is_not_biological_time():
    assert make(date_available=True, elapsed_time_observable=False).qualification != 'CONTINUOUS_TIME_OBSERVABLE'

def test_phase_does_not_create_delta_t():
    q=make(ordered_phase_available=True)
    assert q.elapsed_time_observable is False

def test_unknown_metadata_stays_unknown():
    q=make(data_available=False); assert q.qualification == 'UNKNOWN'

def test_missing_treatment_blocks_continuous_time():
    q=make(elapsed_time_observable=True, treatment_metadata_available=False)
    assert q.qualification == 'LONGITUDINAL_BUT_INSUFFICIENT'

def test_treatment_is_not_time_proxy():
    q=make(elapsed_time_observable=False, ordered_phase_available=False)
    assert temporal_dataset_gate(q)['continuous_time_available'] is False
