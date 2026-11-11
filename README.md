# Territorial Cerebrovascular Atlas

Analysis code and frozen result tables for:

**Arterial territory as a molecular and imaging coordinate of human cerebrovascular disease**

## Contents

- nalysis/ — analysis scripts for molecular Class B contrasts, imaging geography (ATLAS / ISLES / ADNI), MultiVINE genetics, developmental and Circle of Willis tests, and related sensitivity analyses
- 
esults/ — frozen numerical outputs (CSV/JSON/Markdown reports) corresponding to the manuscript analyses

## Requirements

Python 3.10+ recommended. See 
equirements.txt for core dependencies (scanpy/anndata, scikit-learn, statsmodels, and related scientific Python packages).

## Data

Analyses use publicly available and institutionally controlled resources described in the manuscript (Garcia ROSMAP vascular snRNA-seq, Liu arterial-territory atlas, ATLAS, ISLES, ADNI, MultiVINE, and related cohorts). Users should obtain each dataset through its official access channel and place files under a local data/ layout consistent with path helpers in nalysis/.

## Reproducibility

Primary analysis entry points live under:

- nalysis/moonshot/ — genetics, transport sensitivity, developmental fate, CoW prediction, manifold/borderzone, ecological models
- nalysis/aim4/ — territory-normalized lesion and WMH geography
- nalysis/tracks/ — imaging track analyses including ISLES’24 CoW features
- nalysis/aca_cohort_hunt/, nalysis/garcia/, nalysis/yang/ — cohort assembly and supporting contrasts

Locked numeric outputs are under the matching folders in 
esults/.

## Citation

If you use this code or frozen outputs, please cite the accompanying manuscript and the original data sources listed therein.

## License

Analysis code is provided for academic reuse. Upstream datasets remain under their original licenses and data-use terms.
