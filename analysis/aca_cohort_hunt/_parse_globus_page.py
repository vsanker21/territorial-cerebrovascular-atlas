from pathlib import Path
import re
import urllib.request

url = "https://research.libd.org/globus/"
dest = Path(r"F:\Cerebrovascular SCS\data\raw\aca_cohort_hunt\libd_annotations\globus_live.html")
req = urllib.request.Request(url, headers={"User-Agent": "CerebrovascularTerritorialAtlas/1.0"})
dest.write_bytes(urllib.request.urlopen(req, timeout=60).read())
html = dest.read_text(encoding="utf-8", errors="replace")
print("bytes", dest.stat().st_size)
for m in re.finditer(r'href="(https://app\.globus\.org[^"]+)"', html):
    print("LINK", m.group(1))
# also file-manager / collection links
for m in re.finditer(r'https://[^"\s<>]+globus[^"\s<>]*', html, re.I):
    print("URL", m.group(0)[:200])
text = re.sub(r"<[^>]+>", " ", html)
text = re.sub(r"\s+", " ", text)
for key in ["dACC", "spatialdACC", "spatial_dACC", "jhpce#", "UUID"]:
    i = text.lower().find(key.lower())
    if i >= 0:
        print(key, "->", text[max(0, i - 100) : i + 180])
