# kBET Analysis

## Data
Sc_04 dataset: https://doi.org/10.6084/m9.figshare.30480548.v1

**Note:** Update file paths before running.

## Notebooks

### Data Preparation
- `save_celltype_data.ipynb` - Extracts required fields from raw data

### scib kBET Analysis
- `kbet_scib_cell_type_1.ipynb` - Cell type 1 analysis using scib's kBET
  - Uses largest connected component only (smaller component with 120 cells < 3×k₀=210 discarded)
  - Result: kBET = 1.00 (expected)
  
- `kbet_scib_cell_type_2.ipynb` - Cell type 2 analysis using scib's kBET
  - Same connected component approach
  - Result: kBET = 0.26 (unexpected)

### Original kBET Analysis
- `kbet_original_cell_type_1.ipynb` - Cell type 1 with original kBET
  - Uses scib's neighborhood construction (k₀=70) but no connected components
  - Result: kBET = 1.00 (expected)
  
- `kbet_original_cell_type_2.ipynb` - Cell type 2 with original kBET
  - Same approach as cell type 1
  - Result: kBET = 1.00 (expected)

### Validation (Largest Component Only)
- `kbet_original_cell_type_1_largest_comp.ipynb` 
  - Runs original kBET on only the largest connected component of cell type 1
  - Replicates scib's approach of component filtering
  - Result: kBET = 1.00 (matches scib result)
  
- `kbet_original_cell_type_2_largest_comp.ipynb`
  - Runs original kBET on only the largest connected component of cell type 2
  - Replicates scib's approach of component filtering
  - Result: kBET = 0.26 (matches scib result)
 
## Key Findings
Original kBET yields expected scores (1.00) for both cell types, while scib's connected component approach produces unexpected results for cell type 2 (0.26).
