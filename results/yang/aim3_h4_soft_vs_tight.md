# Soft vs tight gate — Aim3 / H4 comparison

## Retention
- Soft: 143,148 -> Tight: 65,382 (45.7%)
- Soft endo/peri: 88,800 / 36,353
- Tight endo/peri: 44,059 / 15,121

## LINGO1 rank in Aim3 DE (lower/absent = cleaner)

| Stratum | Soft gate | Tight gate |
|---|---|---|
| endo/control | rank 2, logFC=3.13, padj=7.89e-179 | not in top list |
| endo/AD | rank 11, logFC=2.35, padj=1.54e-68 | not in top list |
| peri/control | rank 6, logFC=2.36, padj=5.13e-66 | not in top list |
| peri/AD | rank 33, logFC=1.81, padj=1.09e-30 | not in top list |

## H4 AUROC

| Gate | Control CV | AD transfer | Delta |
|---|---:|---:|---:|
| Soft | 0.870 | 0.665 | -0.204 |
| Tight | 0.878 | 0.695 | -0.182 |
