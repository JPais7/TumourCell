"""Conservative malignant-cell observability audit for GSE205472."""
RULES=('published_annotation','cnv_inference','epithelial_marker','immune_stromal_exclusion')
def classify_rule(*, published_annotation=False, cnv_available=False, epithelial=False, outcome_used=False, treatment_specific=False):
    if outcome_used or treatment_specific: return 'UNRESOLVED'
    if published_annotation and cnv_available: return 'MALIGNANT_PROBABLE'
    if cnv_available: return 'MALIGNANT_UNCERTAIN'
    if epithelial: return 'MALIGNANT_UNCERTAIN'
    return 'UNRESOLVED'
def phase_gate(status):
    return status=='MALIGNANT_CELL_OBSERVABLE_VERIFIED'
def sensitivity_status(statuses):
    return {'rules_evaluated':list(statuses),'promote':False if any(s in ('UNRESOLVED','MALIGNANT_UNCERTAIN') for s in statuses) else False}
