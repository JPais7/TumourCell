from src.malignant.gse205472 import classify_rule, phase_gate, sensitivity_status
def test_epithelial_is_not_malignant(): assert classify_rule(epithelial=True)=='MALIGNANT_UNCERTAIN'
def test_tumour_sample_is_not_cell_ground_truth(): assert classify_rule()=='UNRESOLVED'
def test_outcome_cannot_define_malignancy(): assert classify_rule(outcome_used=True,cnv_available=True)=='UNRESOLVED'
def test_treatment_specific_rule_is_leakage(): assert classify_rule(treatment_specific=True,cnv_available=True)=='UNRESOLVED'
def test_cnv_is_not_automatic_ground_truth(): assert classify_rule(cnv_available=True)=='MALIGNANT_UNCERTAIN'
def test_partial_never_unlocks_phase5(): assert phase_gate('MALIGNANT_PROBABLE') is False
def test_sensitivity_does_not_choose_best_effect(): assert sensitivity_status(['MALIGNANT_UNCERTAIN','MALIGNANT_PROBABLE'])['promote'] is False
