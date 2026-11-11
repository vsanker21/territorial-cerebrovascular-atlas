"""Print Garcia SI workbook structure for audit refinement."""
from pathlib import Path
import pandas as pd

path = Path(__file__).resolve().parents[2] / "data" / "raw" / "garcia" / "Garcia_SuppTable1.xlsx"
xl = pd.ExcelFile(path)
print("sheets:", xl.sheet_names)
for s in xl.sheet_names:
    df = pd.read_excel(path, sheet_name=s, header=None)
    print("\n===", s, "shape", df.shape, "===")
    print(df.head(20).to_string())
    print("--- dtypes/cols with header=0 ---")
    df0 = pd.read_excel(path, sheet_name=s)
    print(list(df0.columns))
    print(df0.head(5).to_string())
