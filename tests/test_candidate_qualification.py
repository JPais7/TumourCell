from src.data.candidate_qualification import qualify_candidate, aggregator_label_is_evidence, candidate_rank

def base(**kw):
    d=dict(patient_ids_public=True,repeated_patients_verified=True,malignant_cells_available=True,direct_biological_time=False,discrete_timepoints=False,metadata_access='PUBLIC'); d.update(kw); return d
def test_public_repeated_discrete_candidate(): assert qualify_candidate(base(discrete_timepoints=True))=='DISCRETE_LONGITUDINAL_OBSERVABLE'
def test_continuous_requires_direct_time(): assert qualify_candidate(base(direct_biological_time=True))=='CONTINUOUS_TIME_OBSERVABLE'
def test_public_claim_without_patient_mapping_fails(): assert qualify_candidate(base(patient_ids_public=False))=='NOT_VERIFIED'
def test_missing_malignant_compartment_fails(): assert qualify_candidate(base(malignant_cells_available=False))=='NOT_SUFFICIENT_FOR_MALIGNANT_DYNAMICS'
def test_phase_is_not_discrete_without_verified_mapping(): assert qualify_candidate(base(patient_ids_public=False,discrete_timepoints=True))=='NOT_VERIFIED'
def test_aggregator_label_cannot_override_original(): assert aggregator_label_is_evidence('longitudinal',False) is False
def test_controlled_candidate_not_best(): assert candidate_rank(base(metadata_access='CONTROLLED',patient_ids_public=False))=='PARTIAL'
def test_future_metadata_does_not_change_missing_critical_field(): assert qualify_candidate(base(repeated_patients_verified=False,direct_biological_time=True))=='NOT_VERIFIED'
def test_discrete_verified_pairing_does_not_require_continuous_time(): assert qualify_candidate(base(discrete_timepoints=True,direct_biological_time=False))=='DISCRETE_LONGITUDINAL_OBSERVABLE'
def test_bulk_or_partial_malignant_evidence_cannot_be_primary(): assert qualify_candidate(base(malignant_cells_available='PARTIAL'))=='NOT_SUFFICIENT_FOR_MALIGNANT_DYNAMICS'
