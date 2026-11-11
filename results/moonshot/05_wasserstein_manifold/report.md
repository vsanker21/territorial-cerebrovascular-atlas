# Moonshot 05 — Wasserstein territorial manifold

**Honesty lock:** Theory figure — pair with genetics/imaging; not standalone claim

## Pairwise W2
| t1   | t2   |       W2 |
|:-----|:-----|---------:|
| ACA  | ACA  | 17.2404  |
| ACA  | MCA  | 30.0175  |
| MCA  | MCA  |  2.15671 |

## H2 OT borderzone
```json
{
  "W2_SFG_ACA": 55.71717508981522,
  "W2_SFG_MCA": 9.971599020430665,
  "W2_SFG_geodesic_mid": 27.039355913022305,
  "W2_ACA_MCA": 52.24239868940731,
  "OT_mixture_FVE": 0.5435929430387354,
  "prior_H2_FVE": 0.94,
  "interpretation": "If OT_mixture_FVE high, SFG sits near ACA\u2013MCA geodesic (mixture/borderzone), consistent with H2; SFG is NOT clean ACA."
}
```

## Conformal labels
```json
{
  "accuracy": 0.8166259168704156,
  "set_coverage": 0.921760391198044,
  "mean_set_size": 1.5525672371638142,
  "n_test": 409
}
```