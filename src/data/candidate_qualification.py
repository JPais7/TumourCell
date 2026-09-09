"""Deterministic dataset-candidate qualification for Phase 4.6."""
CRITICAL=("patient_ids_public","repeated_patients_verified","malignant_cells_available")

def qualify_candidate(candidate):
    """Return discrete/continuous status without arbitrary numeric scoring."""
    if candidate.get('patient_ids_public') is not True or candidate.get('repeated_patients_verified') is not True:
        return 'NOT_VERIFIED'
    if candidate.get('malignant_cells_available') is not True:
        return 'NOT_SUFFICIENT_FOR_MALIGNANT_DYNAMICS'
    if candidate.get('direct_biological_time') is True:
        return 'CONTINUOUS_TIME_OBSERVABLE'
    if candidate.get('discrete_timepoints') is True:
        return 'DISCRETE_LONGITUDINAL_OBSERVABLE'
    return 'ORDERED_PHASE_ONLY'

def aggregator_label_is_evidence(label, original_verified):
    return bool(original_verified) and label in ('longitudinal','continuous_time')

def candidate_rank(candidate):
    """Qualitative ranking; critical missing fields cannot be compensated."""
    status=qualify_candidate(candidate)
    if status=='CONTINUOUS_TIME_OBSERVABLE': return 'BEST_CANDIDATE'
    if status=='DISCRETE_LONGITUDINAL_OBSERVABLE': return 'PROMISING'
    if candidate.get('metadata_access')=='CONTROLLED': return 'PARTIAL'
    return 'INSUFFICIENT'
