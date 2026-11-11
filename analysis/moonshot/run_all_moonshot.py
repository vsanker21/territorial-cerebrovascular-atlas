"""Run all moonshot analyses 01–07 and write synthesis."""
from __future__ import annotations

import json
import traceback
from datetime import datetime, timezone
from pathlib import Path

from _paths import HONESTY, RES


def main():
    import importlib

    modules = [
        ("01", "01_gwas_multivine"),
        ("02", "02_enrichment_transport"),
        ("03", "03_developmental_fate"),
        ("04", "04_cow_iv_conformal"),
        ("05", "05_wasserstein_manifold"),
        ("06", "06_bayesian_ecological"),
        ("07", "07_foundation_probes"),
    ]
    summaries = {}
    for key, modname in modules:
        print("\n" + "=" * 60, flush=True)
        print("RUNNING", modname, flush=True)
        try:
            mod = importlib.import_module(modname)
            # reload if re-run
            mod = importlib.reload(mod)
            summaries[key] = mod.run()
            summaries[key]["_status"] = "ok"
        except Exception as e:
            summaries[key] = {
                "_status": "error",
                "error": str(e),
                "traceback": traceback.format_exc(),
                "honesty_lock": HONESTY.get(key),
            }
            print("ERROR", modname, e, flush=True)
            (RES / f"{key}_ERROR.txt").write_text(traceback.format_exc(), encoding="utf-8")

    # Venue gate synthesis
    s1 = summaries.get("01", {})
    s4 = summaries.get("04", {})
    s3 = summaries.get("03", {})
    falsifier = (s1.get("falsifier") or {}) if isinstance(s1, dict) else {}
    genetics_ok = falsifier.get("gwas_beats_null") is True or (
        isinstance(s1.get("top_enrichments"), list) and len(s1.get("top_enrichments") or []) > 0
    )
    trackc = s4.get("track_C_topology_support")
    trackb = s3.get("track_B_language")

    if genetics_ok and trackc == "PASS_EXPLORATORY":
        venue = "CONDITIONAL_TRY_NAT_NEUROSCI — genetics + topology OOS both non-null; still need manuscript-grade packaging + ACC/RADC for full bar"
    elif genetics_ok:
        venue = "BRAIN_OR_NAT_COMMUN — genetics spine present; Track C OOS weak; ACC/RADC still upgrade path"
    else:
        venue = "BRAIN — public moonshot stack complete with honest nulls; do not force Nat Neurosci"

    synthesis = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "venue_gate": venue,
        "honesty_locks": HONESTY,
        "summaries": summaries,
    }
    (RES / "SYNTHESIS.json").write_text(json.dumps(synthesis, indent=2, default=str), encoding="utf-8")

    lines = [
        "# Moonshot stack synthesis — public-data Nat Neurosci path",
        "",
        f"**UTC:** {synthesis['timestamp_utc']}",
        f"**Venue gate:** {venue}",
        "",
        "## Per-analysis status",
    ]
    for key, modname in modules:
        st = summaries.get(key, {})
        lines.append(f"- **{key} ({modname})**: `{st.get('_status')}` — {st.get('honesty_lock', HONESTY.get(key, ''))}")
        if st.get("_status") == "error":
            lines.append(f"  - ERROR: {st.get('error')}")
    lines += [
        "",
        "## Key locks",
        f"- #1 GWAS beats null: `{falsifier.get('gwas_beats_null')}`",
        f"- #3 Track B language: `{trackb}`",
        f"- #4 Track C topology: `{trackc}`",
        f"- #2 is SENSITIVITY_NOT_H1c always",
        f"- #6 is ECOLOGICAL_NOT_CONFIRMATORY_H3 always",
        "",
        "See per-analysis `report.md` under `results/moonshot/0*/`.",
    ]
    (RES / "SYNTHESIS.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nSYNTHESIS written", RES / "SYNTHESIS.md")
    print("VENUE GATE:", venue)


if __name__ == "__main__":
    main()
