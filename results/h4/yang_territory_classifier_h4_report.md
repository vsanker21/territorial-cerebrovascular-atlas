# H4 first-pass — Yang endothelial territory classifier

> **SUPERSEDED for H4 claims.** Formal lock is `yang_h4_formal_status.md` → `H4_NOT_SUPPORTED_ON_YANG_ALONE` (within-donor AUROC stays high in AD; Δ positive). This first-pass cell-CV transfer AUROC drop is **exploratory / not the confirmatory H4 metric**.

Trained on **control** endo: PCA hippocampus vs ACA/MCA borderzone SFG.
Transferred to **AD** endo without retraining.

| Metric | Value |
|---|---:|
| control_cv_auroc | 0.8698 |
| control_cv_accuracy | 0.7909 |
| ad_transfer_auroc | 0.6654 |
| ad_transfer_accuracy | 0.6174 |
| delta_auroc_ad_minus_control | -0.2043 |
| control_mean_correct_score | 0.7448 |
| ad_mean_correct_score | 0.5949 |
| n_control | 16000 |
| n_ad | 16000 |
| n_features | 1501 |
| lineage | endo |
| classes | ACA_MCA_borderzone_SFG vs PCA_hippocampus |
| interpretation | H4 predicts AD transfer AUROC < control CV AUROC (territorial identity degradation). First-pass exploratory. |

**Delta AUROC (AD - control CV) = -0.2043**

Negative delta is consistent with H4 (degradation); treat as exploratory until donor-aware CV and lineage-cleaned markers are locked.
