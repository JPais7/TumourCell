Interpretation: each RNA-seq run is treated as one cell/nucleus library because the ENA sample_alias identifies a patient-specific cell.
Timepoints remain candidate labels: 0cell -> pre_candidate, 2cell -> mid_candidate, OPcell -> unknown.
Clone identity is not inferred by this file; it must be linked from the DNA-cell tables/clone assignments.
The companion file `kim_pilot_cell_clone_linkage.csv` keeps clone identity as
`unresolved` for all 60 cells. It adds the patient-level clonal outcome from
Table S4 (`Extinction` or `Persistence`), which is not a cell-level clone
label and must not be used as one.
