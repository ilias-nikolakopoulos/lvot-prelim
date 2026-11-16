"""Build outputs/email/LVOT_findings_email.eml (HTML with inline figures) and a plain-text copy.
Open the .eml in Apple Mail / Outlook, then forward or copy into a new message. Not sent anywhere.
Usage: python outputs/email/build_email.py
"""
from email.message import EmailMessage
from email.utils import make_msgid
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUBJECT = "LVOT size and shape in UK Biobank average hearts: findings so far"

FIGS = [
    ("fig1.png", "Figure 1. How the LVOT is cut: planes parallel to the aortic annulus, 2 mm and 4 mm below it. "
                 "At 4 mm the outline is open on the mitral side (no leaflets in the meshes) and is closed with a straight chord."),
    ("fig2.png", "Figure 2. Men's annulus and LVOT are about 3 mm wider than women's at every age (same-age, same-BMI pairs)."),
    ("fig3.png", "Figure 3. The annulus is nearly round; the LVOT just below it is clearly oval, more so deeper and more so in women."),
    ("fig4.png", "Figure 4. Ellipse formulas get the LVOT area within a few percent; the echo circle from the PLAX diameter "
                 "underestimates it by about 30%, more in women."),
    ("fig5.png", "Figure 5. Features that go with a more oval LVOT (all values as z-scores for sex, age and BMI). "
                 "Note: 'Aorto-mitral plane angle' here is the angle between the ring normals, i.e. 180° minus the clinical "
                 "aorto-mitral angle, so a negative ρ means a wider clinical angle goes with a more oval LVOT."),
]

PARAS = [
    ("p", "Hi all,"),
    ("p", "Here is a short summary of where the LVOT side project stands. We measured the left ventricular outflow tract "
          "(LVOT) in the 1,423 representative UK Biobank heart meshes released by Ugurlu et al. (PLoS One 2025). "
          "Each mesh is the average end-diastolic heart for one sex, age (46–83) and BMI, so the results describe typical "
          "hearts, not individuals. All 1,423 meshes were processed with no errors."),
    ("h", "How the LVOT is measured"),
    ("p", "The LVOT has no natural lower edge, so every cut is anchored to the aortic annulus (where the aortic valve "
          "meets the LV wall in the mesh). Cuts are flat planes parallel to the annulus, 2 mm below it (main result: a "
          "closed ring of wall) and 4 mm below it (second result). The meshes have no valve leaflets, so below about 3 mm "
          "the mitral side of the LVOT is open and is bridged with a straight line (Figure 1)."),
    ("fig", 0),
    ("h", "1. Size: men's outflow is about 3 mm wider at every age"),
    ("p", "Comparing men and women of the same age and BMI (622 pairs), the aortic annulus and LVOT are 2.6–3.0 mm "
          "(15–16%) wider in men in every age group. Median LVOT diameter at 2 mm is 18.5 mm in women and 21.4 mm in men. "
          "The gap does not narrow with age, and it stays at about 6% after indexing to LV size, so LV size does not fully "
          "explain it (Figure 2)."),
    ("fig", 1),
    ("h", "2. Shape: the LVOT is oval, more so in women"),
    ("p", "The annulus is nearly round (major/minor diameter 1.14 in both sexes). Just below it the LVOT is clearly oval: "
          "1.47 in women vs 1.40 in men at 2 mm, and 1.53 vs 1.43 at 4 mm. Women's LVOT is more elliptical at every age, "
          "and this does not change with age. The short axis runs from the septum to the aorto-mitral side, close to the "
          "echo parasternal long-axis line (Figure 3)."),
    ("fig", 2),
    ("h", "3. The echo circle underestimates LVOT area by about 30%"),
    ("p", "Ellipse estimates from two diameters are within 7% of the true outline area. The standard echo "
          "calculation, a circle from the PLAX diameter, measures the short axis and underestimates the area by about "
          "33% in women and 29% in men at 2 mm. The extra 4–6 percentage points in women would carry into stroke volume and "
          "continuity-equation valve area (Figure 4)."),
    ("fig", 3),
    ("h", "4. What goes with a more oval LVOT: valve-plane geometry, not wall thickness"),
    ("p", "We added regional wall thickness (12 segments plus the subaortic septum) and the aorto-septal angle, and "
          "compared each mesh with others of the same sex, age and BMI. The strongest link is the aorto-mitral angle: "
          "a wider (flatter) clinical aorto-mitral angle goes with a more oval LVOT. Across all 998 meshes with normal "
          "walls, a 1 SD wider angle goes with about 0.4 SD more ellipticity, and this holds after adjusting for annulus "
          "size and wall thickness. A smaller aortic annulus also goes with a more oval LVOT. Wall thickness itself shows "
          "almost no link, and meshes with a normal septum but thicker walls elsewhere are no different (Figure 5)."),
    ("fig", 4),
    ("p", "A quick PubMed review fits this: a wider aorto-mitral angle goes with a larger LV and mitral annulus, while "
          "a narrower angle goes with older age, a small hyperdynamic LV, aortic elongation, and obstructive HCM and SAM. "
          "So if LVOT shape were an early HCM marker, our data suggest a rounder LVOT with a narrow aorto-mitral angle, "
          "not a more oval one. That needs testing in individual scans."),
    ("h", "Main caveats"),
    ("ul", ["End-diastole only: the dataset has no systolic meshes, while clinical LVOT measurements are usually mid-systolic.",
            "No valve leaflets: the LVOT is a closed tube for only 2–5 mm; deeper sections use a straight chord.",
            "Averaged hearts: each mesh averages many people, which shrinks the spread (e.g. aorto-mitral angle SD about "
            "2.6° vs about 10° in patients). P-values are descriptive.",
            "Cross-sectional: age differences are between groups, not changes over time."]),
    ("h", "Possible next steps"),
    ("ul", ["Repeat the key measures (LVOT ellipticity and aorto-mitral angle) in individual UK Biobank CMR scans, in "
            "systole and diastole.",
            "Test whether a narrow aorto-mitral angle with a rounder LVOT is more common in HCM or in carriers of HCM variants."]),
    ("p", "Code, tables and full methods are in the private GitHub repository (ilias-nikolakopoulos/lvot-prelim). Happy to "
          "go through any of it."),
    ("p", "Best wishes,<br>[Name]"),
]


def html():
    cids = [make_msgid(domain="lvot.local")[1:-1] for _ in FIGS]
    out = ['<div style="font-family:Arial,Helvetica,sans-serif;font-size:14px;line-height:1.5;color:#102D40;max-width:760px">']
    for kind, v in PARAS:
        if kind == "p":
            out.append(f"<p>{v}</p>")
        elif kind == "h":
            out.append(f'<p style="margin:20px 0 6px"><b>{v}</b></p>')
        elif kind == "ul":
            out.append("<ul>" + "".join(f"<li>{x}</li>" for x in v) + "</ul>")
        else:
            f, cap = FIGS[v]
            out.append(f'<p><img src="cid:{cids[v]}" alt="{cap}" style="width:100%;max-width:760px;border:0"><br>'
                       f'<span style="font-size:12px;color:#5A6B77">{cap}</span></p>')
    return "\n".join(out) + "</div>", cids


def text():
    out = []
    for kind, v in PARAS:
        if kind in ("p", "h"):
            out.append(v.replace("<br>", "\n"))
        elif kind == "ul":
            out.append("\n".join(f"- {x}" for x in v))
        else:
            out.append(f"[{FIGS[v][1]}  (attached: {FIGS[v][0]})]")
    return "\n\n".join(out) + "\n"


body, cids = html()
msg = EmailMessage()
msg["Subject"] = SUBJECT
msg["To"] = ""
msg["X-Unsent"] = "1"  # Outlook / Apple Mail open it as an editable draft
msg.set_content(text())
msg.add_alternative(body, subtype="html")
html_part = msg.get_payload()[1]
for (f, _), cid in zip(FIGS, cids):
    html_part.add_related((HERE / "img" / f).read_bytes(), "image", "png", cid=f"<{cid}>", filename=f)
(HERE / "LVOT_findings_email.eml").write_bytes(bytes(msg))
(HERE / "LVOT_findings_email.txt").write_text(f"Subject: {SUBJECT}\n\n" + text())
preview = body
for (f, _), cid in zip(FIGS, cids):
    preview = preview.replace(f"cid:{cid}", f"img/{f}")
(HERE / "LVOT_findings_email.html").write_text(
    "<!doctype html><meta charset='utf-8'><title>LVOT findings email</title><body style='background:#fff'>"
    + preview + "</body>")
print("wrote LVOT_findings_email.eml, .txt, .html")
