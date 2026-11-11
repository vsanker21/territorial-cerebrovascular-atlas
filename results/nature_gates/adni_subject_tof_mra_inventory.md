# ADNI subject-level TOF/MRA inventory — Track C gate

**Date:** 2026-08-18
**Verdict:** `ADNI_SUBJECT_TOF_MRA__NO_GO_ON_DISK_AND_PROTOCOL`

## Finding
- MRI key rows scanned: **92,246** (3,195 subjects).
- Text hits for MRA/TOF/angio/TOF synonyms: **0**.
- Present `series_type` values include T1W, FLAIR, ASL, dMRI, BOLD, etc. — **no MRA/TOF class**.
- Aim 4 WMH nest n=155: MRA columns are **atlas-constant**, not per-subject angio.

## Protocol context
ADNI core MRI (ADNI1–4 documentation) emphasizes T1, FLAIR, T2*, dMRI, ASL, tf-fMRI — **not** TOF-MRA. Local CSV search across all ADNI tables found no TOF/MRA modality fields (only `MRACQUISITIONTYPEFLAIR`, unrelated).

## Scientific implication (Nature Track C)
Track C subject-level architecture mediation cannot be completed on ADNI WMH. Keep ADNI for territory-normalized WMH geography only. Pursue stroke/aneurysm cohorts with subject TOF+lesions, or OpenNeuro TOF methods pilot.

## Do not claim
- That atlas-constant MRA null proves or refutes subject-level architecture mediation.
- That ADNI WMH geography alone is a complete Nature mechanism package.

## Next
- OpenNeuro ds005096 TOF-MRA (methods/calibration; no stroke lesions)
- Stroke center cohorts with TOF/CTA + MNI lesions (request/collaborate)
- Do not re-interpret atlas MRA null as subject architecture
