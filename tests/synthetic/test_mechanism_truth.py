from src.clones.clone_state import assess_observed_shift,MechanismStatus
def test_pure_selection_requires_all_evidence():
 a=assess_observed_shift(clone_evidence=True,assignment_quality=.95,temporal_overlap=True,abundance_sufficient=True,sampling_comparable=True,clone_frequency_changed=True)
 assert a.selection==MechanismStatus.SUPPORTED.value; assert a.plasticity!=MechanismStatus.SUPPORTED.value
def test_pure_plasticity_requires_within_clone_state_change():
 a=assess_observed_shift(clone_evidence=True,assignment_quality=.95,temporal_overlap=True,sampling_comparable=True,within_clone_state_change=True)
 assert a.plasticity==MechanismStatus.SUPPORTED.value; assert a.selection!=MechanismStatus.SUPPORTED.value
def test_mixed_is_not_forced_to_one_mechanism():
 a=assess_observed_shift(clone_evidence=True,assignment_quality=.95,temporal_overlap=True,abundance_sufficient=True,sampling_comparable=True,clone_frequency_changed=True,within_clone_state_change=True)
 assert a.selection==a.plasticity==MechanismStatus.SUPPORTED.value; assert a.identifiability_status==MechanismStatus.INDETERMINATE.value
