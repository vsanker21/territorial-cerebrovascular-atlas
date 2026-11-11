# H4 first-pass — Yang endothelial territory classifier (gate v2 tight)

Trained on **control** endo (tight): PCA hippocampus vs ACA/MCA borderzone SFG.
Transferred to **AD** endo without retraining.

| Metric | Value |
|---|---:|
| control_cv_auroc | 0.8777 |
| control_cv_accuracy | 0.8028 |
| ad_transfer_auroc | 0.6952 |
| ad_transfer_accuracy | 0.6574 |
| delta_auroc_ad_minus_control | -0.1825 |
| control_mean_correct_score | 0.7541 |
| ad_mean_correct_score | 0.6244 |
| n_control | 16000 |
| n_ad | 11353 |
| n_features | 1500 |
| lineage | endo |
| classes | ACA_MCA_borderzone_SFG vs PCA_hippocampus |
| interpretation | H4 predicts AD transfer AUROC < control CV AUROC (territorial identity degradation). First-pass exploratory. |

**Delta AUROC (AD - control CV) = -0.1825**

Negative delta is consistent with H4 (degradation); still exploratory until donor-aware CV.
