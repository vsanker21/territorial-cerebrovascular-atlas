# Winkler UCSC — zonation vs territory

**Verdict:** `ZONATION_PROGRAMS_PRESENT_TERRITORY_NO_GO`

- A–V-like EC clusters present (Art/Cap/Venous/Venule).
- **No region / TerritoryVoxels metadata** → cannot test residual territory signal within Winkler.
- Use as zonation covariate library / AVM disease exploratory only.
- **Caveat:** near-perfect Cap vs Art and AVM vs CTRL cell CV AUROCs are expected when labels are discrete clusters / sample-confounded strata — treat as **label sanity**, not biological proof of H4.

CTRL Cap vs Arterial cell CV AUROC: **1.000**
AVM vs CTRL (am-endothelial) cell CV AUROC: **1.000** (exploratory H4-adjacent, not territory).
