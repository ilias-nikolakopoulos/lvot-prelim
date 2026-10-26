"""STEP 7f (part 2) — Add the on-click build to outputs/slide_ellipticity_animated.pptx.

pptxgenjs cannot write animations, so this adds a <p:timing> block to slide 1: on the first click,
"Headline sex" and "Chart build 2" appear and "Chart build 1" disappears (all at once).
Shapes are found by the objectName set in scripts/20_ellipticity_slide.js. Run after that script.
Usage: python scripts/20_add_click_build.py
"""
import re
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.provenance import log_run  # noqa: E402

PPTX = ROOT / "outputs/slide_ellipticity_animated.pptx"
SLIDE = "ppt/slides/slide1.xml"


def effect(ctn, spid, cls, node, visible):
    """One 'Appear' (entr) or 'Disappear' (exit) effect on shape spid."""
    return (f'<p:par><p:cTn id="{ctn}" presetID="1" presetClass="{cls}" presetSubtype="0" fill="hold" '
            f'nodeType="{node}"><p:stCondLst><p:cond delay="0"/></p:stCondLst><p:childTnLst><p:set><p:cBhvr>'
            f'<p:cTn id="{ctn + 1}" dur="1" fill="hold"><p:stCondLst><p:cond delay="0"/></p:stCondLst></p:cTn>'
            f'<p:tgtEl><p:spTgt spid="{spid}"/></p:tgtEl><p:attrNameLst><p:attrName>style.visibility</p:attrName>'
            f'</p:attrNameLst></p:cBhvr><p:to><p:strVal val="{"visible" if visible else "hidden"}"/></p:to>'
            f'</p:set></p:childTnLst></p:cTn></p:par>')


with zipfile.ZipFile(PPTX) as z:
    parts = {n: z.read(n) for n in z.namelist()}
xml = parts[SLIDE].decode("utf8")
ids = {name: int(i) for i, name in re.findall(r'<p:cNvPr id="(\d+)" name="([^"]+)"', xml)}
need = ["Headline sex", "Chart build 2", "Chart build 1"]
missing = [n for n in need if n not in ids]
if missing:
    sys.exit(f"shapes not found in {SLIDE}: {missing}; found {sorted(ids)}")
if "<p:timing>" in xml:
    sys.exit("slide already has a <p:timing> block; rebuild the deck first")

clicks = (effect(5, ids["Headline sex"], "entr", "clickEffect", True)
          + effect(7, ids["Chart build 2"], "entr", "withEffect", True)
          + effect(9, ids["Chart build 1"], "exit", "withEffect", False))
timing = ('<p:timing><p:tnLst><p:par><p:cTn id="1" dur="indefinite" restart="never" nodeType="tmRoot"><p:childTnLst>'
          '<p:seq concurrent="1" nextAc="seek"><p:cTn id="2" dur="indefinite" nodeType="mainSeq"><p:childTnLst>'
          '<p:par><p:cTn id="3" fill="hold"><p:stCondLst><p:cond delay="indefinite"/></p:stCondLst><p:childTnLst>'
          '<p:par><p:cTn id="4" fill="hold"><p:stCondLst><p:cond delay="0"/></p:stCondLst><p:childTnLst>'
          + clicks +
          '</p:childTnLst></p:cTn></p:par></p:childTnLst></p:cTn></p:par></p:childTnLst></p:cTn>'
          '<p:prevCondLst><p:cond evt="onPrev" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond></p:prevCondLst>'
          '<p:nextCondLst><p:cond evt="onNext" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond></p:nextCondLst>'
          '</p:seq></p:childTnLst></p:cTn></p:par></p:tnLst></p:timing>')
# <p:timing> goes after cSld / clrMapOvr / transition, before extLst
pos = max(m.end() for m in re.finditer(r"</p:clrMapOvr>|</p:cSld>|</p:transition>", xml))
xml = xml[:pos] + timing + xml[pos:]
parts[SLIDE] = xml.encode("utf8")

tmp = Path(tempfile.mkstemp(suffix=".pptx")[1])
with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
    for n, data in parts.items():  # keeps [Content_Types].xml first, as written by pptxgenjs
        z.writestr(n, data)
shutil.move(tmp, PPTX)
print(f"added click build to {PPTX.relative_to(ROOT)}: shape ids {[ids[n] for n in need]}")
log_run(__file__, output=str(PPTX.relative_to(ROOT)))
