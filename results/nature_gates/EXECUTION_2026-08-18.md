# Nature gates execution — 2026-08-18

Executed all five Nature-path items at highest scientific accuracy. Agent cannot send email or deposit OSF without your credentials.

## 1. OSF prereg v2 (do not overwrite v1.2)

| Item | Path / status |
|---|---|
| Deposit package | `docs/osf_deposit_v2/` |
| Locked copy of v1.2 | `01_prereg_v1.2_LOCKED.md` |
| Amendment v2 | `02_prereg_Nature_amendment_v2.md` |
| Upload script | `deposit_osf.py` (refuses overwrite of same remote name) |
| Blocker | `OSF_TOKEN` unset — script exit 2 |

**Your next step (you selected token deposit):** in PowerShell:

```powershell
$env:OSF_TOKEN = "<token from https://osf.io/settings/tokens>"
$env:OSF_PROJECT = "<5-char project GUID>"
python "F:\Cerebrovascular SCS\docs\osf_deposit_v2\deposit_osf.py"
```

Then ask the agent to re-run if preferred. Manual alternative: `README_DEPOSIT.md`.

## 2. RADC RRA-539487 (phenotype linkage)

| Item | Path |
|---|---|
| Push emails (Stanford + Greg) | `docs/outreach/radc_stanford_dua/email_push_RRA539487_2026-08-18.md` |
| Aim 3 scaffold | Still `GATED_AWAITING_RADC_DUA_DELIVERY` |
| Drop path when delivered | `data/raw/radc/` |

**Send status:** drafts only (per your choice). Scientific role: Tier-2 within-donor phenotype — not identity proof alone.

## 3. Vessel-enriched ACC (only confirmatory DE path)

| Item | Path / verdict |
|---|---|
| Final Wälchli nudge | `docs/outreach/email_Walchli_final_nudge_2026-08-18.md` |
| Timeout | **2026-08-25** |
| GEO / literature hunt | **No public vessel-enriched ACC** |
| Unenriched sgACC cortical snRNA | **NO-GO** for H1c (enrichment mismatch) |
| Yang SFG | Still **NO-GO** as clean ACA |

## 4. Mechanism: subject-level angio/MRA **or** snATAC

| Item | Path / verdict |
|---|---|
| Spec + decision tree | `results/nature_gates/MECHANISM_angio_or_snATAC.md` |
| Track C | Subject-level angio (ADNI TOF inventory next); OpenNeuro aneurysm TOF = methods pilot only |
| Track B | snATAC/MultiVINE only if chromatin works; HOX RNA remains NO-GO |
| Yang MultiVINE ask | `docs/outreach/email_Yang_MultiVINE_snATAC_ask_2026-08-18.md` (confirm corresponding email before send) |

## 5. TRCI/BCI out of Nature package

| Item | Path |
|---|---|
| Explicit exclusion lock | `docs/14_Nature_package_TRCI_BCI_exclusion.md` |
| Checklist | BCI marked OUT of Nature package |

## Scientific claim hygiene (unchanged)

- Confirmatory DE only after enrichment-matched ACA vs MCA.  
- Atlas-constant MRA ≠ subject-level architecture.  
- Endo separability ≠ crest-lineage identity.  
- Track B not lead until chromatin + ACC gates open.
