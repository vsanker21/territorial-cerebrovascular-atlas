from pathlib import Path
import pandas as pd

esm = Path(__file__).resolve().parents[2] / "data" / "raw" / "garcia" / "esm"
for path in sorted(esm.glob("*.xlsx")):
    print("\n########", path.name, path.stat().st_size, "########")
    xl = pd.ExcelFile(path)
    print("sheets:", xl.sheet_names)
    for s in xl.sheet_names[:6]:
        df = pd.read_excel(path, sheet_name=s)
        print(f"  sheet={s!r} shape={df.shape} cols={list(df.columns)[:20]}")
        print(df.head(3).to_string())
