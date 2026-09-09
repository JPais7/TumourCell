from src.validation.real_cohort import temporal_transitions

def obs(patient, phase, treatment=None, row=0):
    treatment = treatment or phase
    return {'row_index': row, 'patient_id': patient, 'sample_id': f'{patient}:observation_{row:03d}',
            'timepoint_id': phase, 'time_value': phase, 'treatment': treatment}

def test_ordered_phase_is_not_equal_time():
    tr, _ = temporal_transitions([obs('p', 0, 'PD1'), obs('p', 2, 'RTPD1', 1)])
    assert tr[0]['delta_t'] is None
    assert tr[0]['delta_t_status'] == 'NOT_AVAILABLE'

def test_irregular_numeric_values_are_preserved_as_ordered_values():
    a, b = obs('p', 7, 'PD1'), obs('p', 21, 'RTPD1', 1)
    tr, _ = temporal_transitions([a, b])
    assert tr[0]['history_time_value'] == 7 and tr[0]['target_time_value'] == 21

def test_duplicate_timestamp_is_rejected():
    tr, failures = temporal_transitions([obs('p', 0), obs('p', 0, row=1), obs('p', 1, row=2)])
    assert tr == []
    assert failures[0]['status'] == 'DUPLICATE_TIMEPOINT'

def test_only_adjacent_forward_transitions_are_created():
    tr, _ = temporal_transitions([obs('p', 0), obs('p', 1, row=1), obs('p', 2, row=2)])
    assert [(x['history_time_value'], x['target_time_value']) for x in tr] == [(0, 1), (1, 2)]

def test_future_response_label_not_in_transition_fields():
    rows = [obs('p', 0), obs('p', 1, row=1)]
    rows[1]['response_label'] = 'R1'
    tr, _ = temporal_transitions(rows)
    assert 'response_label' not in tr[0]

def test_no_future_target_before_history():
    tr, failures = temporal_transitions([obs('p', 2), obs('p', 1, row=1)])
    assert [(x['history_time_value'], x['target_time_value']) for x in tr] == [(1, 2)]
    assert all(x['target_time_value'] > x['history_time_value'] for x in tr)
