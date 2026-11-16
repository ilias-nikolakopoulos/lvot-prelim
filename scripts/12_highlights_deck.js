// STEP 7b — Highlights deck (PowerPoint) built from the current tables and figures.
// Usage: NODE_PATH=<dir with pptxgenjs> node scripts/12_highlights_deck.js
// Output: outputs/LVOT_highlights.pptx   (theme colours applied afterwards with the pptx skill's apply_theme.js)
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");

const ROOT = path.resolve(__dirname, "..");
const T = (f) => path.join(ROOT, "outputs/tables", f);
const FIG = (f) => path.join(ROOT, "outputs/figures", f);

function readCsv(file) {
  const [head, ...rows] = fs.readFileSync(file, "utf8").trim().split("\n");
  const cols = head.split(",");
  return rows.map((r) => Object.fromEntries(r.split(",").map((v, i) => [cols[i], isNaN(+v) || v === "" ? v : +v])));
}
const blocks = readCsv(T("sexdiff_by_age_block.csv"));
const area = readCsv(T("area_methods_summary.csv"));
const ell = readCsv(T("ellipticity_methods.csv"));
const blk = (m) => blocks.filter((r) => r.measure === m);
const AGE = blk("ao_d_area").map((r) => r.age_block.replace("-", "–"));

const THEME = {
  name: "LVOT UKB",
  headFontFace: "Arial",
  bodyFontFace: "Arial",
  colors: {
    dk1: "1E2A33", lt1: "FFFFFF", dk2: "102D40", lt2: "EEF3F6",
    accent1: "B94F00", accent2: "009E73", accent3: "E6C229", accent4: "7B4FA6", accent5: "0072B2", accent6: "D55E00",
    hlink: "0072B2", folHlink: "7B4FA6",
  },
};
const HEX = { women: "D55E00", men: "0072B2", grid: "D9E1E7", axis: "4A5A66", navy: "102D40" };

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.333 x 7.5 in
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "LVOT size and ellipticity in UK Biobank hearts";
pres.author = "Ilias Nikolakopoulos";
const C = pres.SchemeColor;

pres.defineSlideMaster({
  title: "Title dark",
  background: { color: C.text2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.2, w: 11.7, h: 1.5, fontSize: 40, bold: true, color: C.background1, valign: "bottom", align: "left" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 3.85, w: 11.7, h: 1.2, fontSize: 20, color: "CFE0EA", valign: "top", align: "left" }, text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "Title only",
  background: { color: C.background1 },
  margin: [0.5, 0.6, 0.6, 0.6],
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 0.35, w: 12.1, h: 0.9, fontSize: 28, bold: true, color: C.text2, valign: "middle", align: "left" }, text: "" } },
    { text: { text: "LVOT size and ellipticity · UK Biobank representative meshes", options: { x: 0.6, y: 7.0, w: 8, h: 0.3, fontSize: 10, color: "6B7B87" } } },
  ],
  slideNumber: { x: 12.2, y: 7.0, w: 0.5, h: 0.3, fontSize: 10, color: "6B7B87", align: "right" },
});

const chartText = { catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt", legendFontFace: "+mn-lt", titleFontFace: "+mn-lt", dataLabelFontFace: "+mn-lt" };
const chartFrame = () => ({
  ...chartText,
  catAxisLabelColor: HEX.axis, valAxisLabelColor: HEX.axis, catAxisLabelFontSize: 12, valAxisLabelFontSize: 12,
  valGridLine: { color: HEX.grid, size: 0.75 }, catGridLine: { style: "none" },
  showLegend: true, legendPos: "b", legendFontSize: 12, legendColor: HEX.axis,
  showTitle: true, titleFontSize: 15, titleColor: HEX.navy,
});

function card(slide, x, y, w, h, name) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12, fill: { color: C.background2 }, line: { color: C.background2 }, objectName: name });
}
function numberDot(slide, x, y, n, color, name) {
  slide.addShape(pres.shapes.OVAL, { x, y, w: 0.5, h: 0.5, fill: { color }, line: { color }, objectName: name });
  slide.addText(String(n), { x, y, w: 0.5, h: 0.5, fontSize: 16, bold: true, color: C.background1, align: "center", valign: "middle", margin: 0, isTextBox: true, objectName: name + " number" });
}

// ---------------- 1 Title ----------------
pres.addSection({ title: "Overview" });
let s = pres.addSlide({ masterName: "Title dark", sectionTitle: "Overview" });
s.addText("LVOT size and ellipticity in UK Biobank hearts", { placeholder: "title" });
s.addText("Highlights from 1,423 representative end-diastolic meshes, by sex and age · October 2026", { placeholder: "body" });
s.addNotes("Personal exploratory project. Meshes from Ugurlu et al., PLoS One 2025 (Zenodo 10.5281/zenodo.15649643, CC0). All measurements are end-diastolic.");

// ---------------- 2 Key messages ----------------
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Overview" });
s.addText("Four things the meshes show", { placeholder: "title" });
const KM = [
  ["~3 mm", C.accent5, "Men's annulus and LVOT are larger", "2.7–2.9 mm (15–16%) at every age, with age and BMI held fixed; ~6% after indexing to LV size"],
  ["1.14 vs 1.4–1.5", C.accent2, "Round annulus, oval LVOT just below it", "Ellipticity (major/minor) jumps within 2–4 mm of the annulus"],
  ["+0.07 to +0.10", C.accent4, "Women's LVOT is more oval than men's", "At every age, with all three ellipse methods"],
  ["−30%", C.accent1, "Echo circular area underestimates the LVOT", "−33% in women vs −29% in men at 2 mm; any ellipse model is within ~6%"],
];
KM.forEach(([big, col, head, sub], i) => {
  const x = 0.6 + (i % 2) * 6.15, y = 1.55 + Math.floor(i / 2) * 2.65;
  card(s, x, y, 5.95, 2.4, `Key message ${i + 1} card`);
  s.addText(big, { x: x + 0.35, y: y + 0.25, w: 5.3, h: 0.85, fontSize: 40, bold: true, color: col, margin: 0, isTextBox: true, objectName: `Key message ${i + 1} number` });
  s.addText(head, { x: x + 0.35, y: y + 1.1, w: 5.3, h: 0.55, fontSize: 18, bold: true, color: C.text2, margin: 0, isTextBox: true, objectName: `Key message ${i + 1} headline` });
  s.addText(sub, { x: x + 0.35, y: y + 1.65, w: 5.3, h: 0.6, fontSize: 14, color: C.text1, margin: 0, valign: "top", isTextBox: true, objectName: `Key message ${i + 1} detail` });
});
s.addNotes("All values are from the end-diastolic meshes; BMI held fixed by comparing men and women with the same age and BMI (622 pairs).");

// ---------------- 3 Data and approach ----------------
pres.addSection({ title: "Methods" });
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Methods" });
s.addText("1,423 averaged hearts, measured directly on the mesh", { placeholder: "title" });
s.addText([
  { text: "Data", options: { bold: true, color: C.text2, fontSize: 18, breakLine: true } },
  { text: "Representative biventricular meshes, UK Biobank (Ugurlu et al., PLoS One 2025)", options: { bullet: true, breakLine: true } },
  { text: "Each mesh averages participants of the same sex, age (46–83) and BMI", options: { bullet: true, breakLine: true } },
  { text: "761 women, 662 men; 622 same-age, same-BMI pairs hold BMI fixed", options: { bullet: true, breakLine: true } },
  { text: "End-diastole only; valves are flat caps, no leaflets", options: { bullet: true, breakLine: true } },
  { text: " ", options: { breakLine: true } },
  { text: "Output", options: { bold: true, color: C.text2, fontSize: 18, breakLine: true } },
  { text: "Per mesh: annulus, LVOT sections at 0–15 mm, mitral ring, angles, LV volume and length", options: { bullet: true, breakLine: true } },
  { text: "All 1,423 measured, 0 errors; fully scripted and versioned", options: { bullet: true } },
], { x: 0.6, y: 1.5, w: 6.2, h: 5.2, fontSize: 15, color: C.text1, valign: "top", paraSpaceAfter: 6, isTextBox: true, objectName: "Data and output text" });
const STEPS = [
  ["Mesh", "~2 M tetrahedra tagged LV, RV and valve caps"],
  ["LV cavity", "Closed surface of LV wall + mitral and aortic caps"],
  ["Rims and axes", "Aortic and mitral rims, plane fits, LV long axis"],
  ["Plane cuts", "Planes ⟂ aortic axis every 1 mm below the annulus"],
  ["Outline → measures", "Area, diameters, ellipticity, echo-style diameter"],
];
STEPS.forEach(([h, d], i) => {
  const y = 1.55 + i * 1.03;
  card(s, 7.3, y, 5.4, 0.85, `Pipeline step ${i + 1} card`);
  numberDot(s, 7.5, y + 0.17, i + 1, C.text2, `Pipeline step ${i + 1} dot`);
  s.addText([{ text: h, options: { bold: true, color: C.text2, breakLine: true } }, { text: d, options: { color: C.text1, fontSize: 13 } }],
    { x: 8.15, y: y + 0.05, w: 4.45, h: 0.75, fontSize: 15, margin: 0, valign: "middle", isTextBox: true, objectName: `Pipeline step ${i + 1} text` });
});
s.addNotes("Geometry only, no ventricular coordinates. The LV cavity surface splits cleanly into three closed surfaces (outer, RV, LV) because the valve caps close every opening.");

// ---------------- 4 How an outline is measured ----------------
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Methods" });
s.addText("Each LVOT outline is an exact plane cut of the LV wall", { placeholder: "title" });
s.addImage({ path: FIG("qc/outline_steps_female_age_65_bmi_26.png"), x: 0.6, y: 1.35, w: 5.6, h: 5.6, objectName: "Outline steps figure" });
s.addText([
  { text: "From mesh to outline", options: { bold: true, color: C.text2, fontSize: 18, breakLine: true } },
  { text: "A plane perpendicular to the aortic axis cuts the triangulated LV wall", options: { bullet: true, breakLine: true } },
  { text: "Each crossed triangle gives one segment; segments join into a chain", options: { bullet: true, breakLine: true } },
  { text: "Closed chain = closed tube; open chain = the plane hit the mitral opening", options: { bullet: true, breakLine: true } },
  { text: " ", options: { breakLine: true } },
  { text: "From outline to numbers", options: { bold: true, color: C.text2, fontSize: 18, breakLine: true } },
  { text: "Area from the shoelace formula, then diameter = 2√(A/π)", options: { bullet: true, breakLine: true } },
  { text: "Major/minor axes from the moment-equivalent ellipse", options: { bullet: true, breakLine: true } },
  { text: "Only assumption: a straight chord across the mitral-side gap at 4 mm", options: { bullet: true } },
], { x: 6.6, y: 1.5, w: 6.1, h: 5.3, fontSize: 15, color: C.text1, valign: "top", paraSpaceAfter: 6, isTextBox: true, objectName: "Outline method text" });
s.addNotes("Example mesh: woman, 65, BMI 26. Panels: tagged mesh, three closed surfaces, LV endocardium and aortic rim plane, cutting planes, plane-triangle segments, closed chain at 2 mm, open chain at 4 mm, chord closure, shoelace area.");

// ---------------- 5 Short closed tube ----------------
pres.addSection({ title: "Results" });
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Results" });
s.addText("Below the annulus, the LVOT is a closed tube for only 2–5 mm", { placeholder: "title" });
s.addImage({ path: FIG("qc/outline_demo_female_age_65_bmi_26.png"), x: 0.97, y: 1.3, w: 11.4, h: 4.69, objectName: "Outline demo figure" });
const TUBE = [["72%", "of meshes closed for 3 mm"], ["2 mm", "primary level, closed in 99%"], ["4 mm", "secondary level, chord-closed"]];
TUBE.forEach(([b, l], i) => {
  const x = 0.6 + i * 4.1;
  s.addText([{ text: b + "  ", options: { bold: true, fontSize: 22, color: C.accent1 } }, { text: l, options: { fontSize: 14, color: C.text1 } }],
    { x, y: 6.15, w: 3.9, h: 0.5, margin: 0, valign: "middle", isTextBox: true, objectName: `Tube stat ${i + 1}` });
});
s.addNotes("Closed-tube length: 1 mm in 10 meshes, 2 mm in 212, 3 mm in 1,025, 4 mm in 163, 5 mm in 13. Deeper, the section opens towards the mitral valve: the meshes have no anterior mitral leaflet. 6 mm was rejected because large gaps were more common in women.");

// ---------------- 6 Size by sex and age ----------------
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Results" });
s.addText("Men's outflow is ~3 mm wider at every age", { placeholder: "title" });
const lineOpts = (title, lo, hi) => ({
  ...chartFrame(), title, chartColors: [HEX.women, HEX.men], lineSize: 3, lineDataSymbol: "circle", lineDataSymbolSize: 8,
  valAxisMinVal: lo, valAxisMaxVal: hi, valAxisMajorUnit: 1, valAxisTitle: "Diameter, mm", showValAxisTitle: true,
  valAxisTitleFontSize: 12, valAxisTitleColor: HEX.axis, valAxisTitleFontFace: "+mn-lt", catAxisTitle: "Age, years", showCatAxisTitle: true,
  catAxisTitleFontSize: 12, catAxisTitleColor: HEX.axis, catAxisTitleFontFace: "+mn-lt", showValue: false,
});
const series = (m) => [{ name: "Women", labels: AGE, values: blk(m).map((r) => r.mean_F) }, { name: "Men", labels: AGE, values: blk(m).map((r) => r.mean_M) }];
s.addChart(pres.charts.LINE, series("ao_d_area"), { x: 0.6, y: 1.4, w: 4.3, h: 4.6, ...lineOpts("Aortic annulus", 17, 22), objectName: "Annulus diameter chart" });
s.addChart(pres.charts.LINE, series("lvot2_d_area"), { x: 5.0, y: 1.4, w: 4.3, h: 4.6, ...lineOpts("LVOT 2 mm from annulus", 17, 22), objectName: "LVOT diameter chart" });
card(s, 9.6, 1.5, 3.1, 4.4, "Size takeaways card");
s.addText([
  { text: "+2.7 to 2.9 mm", options: { bold: true, fontSize: 24, color: C.accent5, breakLine: true } },
  { text: "men vs women (15–16%), same age and BMI", options: { fontSize: 14, breakLine: true } },
  { text: " ", options: { fontSize: 8, breakLine: true } },
  { text: "Annulus gap narrows slightly with age (2.9 → 2.6 mm)", options: { bullet: true, fontSize: 14, breakLine: true } },
  { text: "LVOT gap stays constant", options: { bullet: true, fontSize: 14, breakLine: true } },
  { text: "Both sexes lose ~0.7 mm from their 50s to late 70s", options: { bullet: true, fontSize: 14 } },
], { x: 9.85, y: 1.7, w: 2.65, h: 4.0, color: C.text1, valign: "top", paraSpaceAfter: 6, isTextBox: true, objectName: "Size takeaways text" });
s.addText("Area-derived diameter, mean of 622 same-age, same-BMI pairs per age block", { x: 0.6, y: 6.1, w: 8.7, h: 0.35, fontSize: 11, color: "6B7B87", isTextBox: true, objectName: "Size chart caption" });
s.addNotes("Paired men − women by age block: annulus gap 2.9 → 2.6 mm (narrows); LVOT 2 and 4 mm gap 2.8–2.9 mm in every block.");

// ---------------- 7 Indexed size ----------------
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Results" });
s.addText("LV size explains most, but not all, of the sex gap", { placeholder: "title" });
s.addChart(pres.charts.LINE, [
  { name: "Annulus, absolute", labels: AGE, values: blk("ao_d_area").map((r) => r.pct_diff) },
  { name: "Annulus ÷ LV volume^⅓", labels: AGE, values: blk("ao_d_area_per_edv13").map((r) => r.pct_diff) },
  { name: "LVOT 2 mm, absolute", labels: AGE, values: blk("lvot2_d_area").map((r) => r.pct_diff) },
  { name: "LVOT 2 mm ÷ LV volume^⅓", labels: AGE, values: blk("lvot2_d_area_per_edv13").map((r) => r.pct_diff) },
], { x: 0.6, y: 1.4, w: 8.2, h: 5.4, ...chartFrame(), title: "Men minus women, % of women's value",
  chartColors: ["B8960F", "B8960F", "009E73", "009E73"], lineDash: ["solid", "dash", "solid", "dash"], lineSize: 3,
  lineDataSymbol: "circle", lineDataSymbolSize: 7, valAxisMinVal: 0, valAxisMaxVal: 18, valAxisMajorUnit: 3,
  showValue: false, valAxisLabelFormatCode: '0"%"', objectName: "Indexed sex gap chart" });
card(s, 9.2, 1.5, 3.5, 5.2, "Indexing takeaways card");
s.addText([
  { text: "~15% → ~6%", options: { bold: true, fontSize: 28, color: C.accent2, breakLine: true } },
  { text: "sex gap after dividing by LV volume^⅓", options: { fontSize: 14, breakLine: true } },
  { text: " ", options: { fontSize: 8, breakLine: true } },
  { text: "Men's outflow stays ~6% larger relative to their LV", options: { bullet: true, fontSize: 14, breakLine: true } },
  { text: "The LV shrinks faster with age than the outflow, so indexed size rises with age, more in women", options: { bullet: true, fontSize: 14, breakLine: true } },
  { text: "Indexed annulus gap narrows: 7.2% → 5.4%", options: { bullet: true, fontSize: 14 } },
], { x: 9.45, y: 1.7, w: 3.05, h: 4.8, color: C.text1, valign: "top", paraSpaceAfter: 6, isTextBox: true, objectName: "Indexing takeaways text" });
s.addNotes("LV cavity volume is the closed LV cavity surface volume (end-diastolic volume proxy): median 116 mL women, 150 mL men.");

// ---------------- 8 Ellipticity ----------------
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Results" });
s.addText("The annulus is round; the LVOT is oval, more so in women", { placeholder: "title" });
s.addImage({ path: FIG("analysis/fig_ellipticity_methods.png"), x: 0.87, y: 1.25, w: 11.6, h: 4.35, objectName: "Ellipticity figure" });
const ELL = [
  ["1.14", "annulus, women and men", C.accent3],
  ["1.46 vs 1.39", "LVOT 2 mm, women vs men", C.accent2],
  ["1.54 vs 1.44", "LVOT 4 mm, women vs men", C.accent4],
  ["±0.03", "agreement between the three ellipse methods", C.text2],
];
ELL.forEach(([b, l, col], i) => {
  const x = 0.6 + i * 3.1;
  card(s, x, 5.8, 2.9, 0.95, `Ellipticity stat ${i + 1} card`);
  s.addText([{ text: b, options: { bold: true, fontSize: 20, color: col, breakLine: true } }, { text: l, options: { fontSize: 12, color: C.text1 } }],
    { x: x + 0.2, y: 5.85, w: 2.6, h: 0.85, margin: 0, valign: "middle", isTextBox: true, objectName: `Ellipticity stat ${i + 1} text` });
});
s.addNotes("Ellipticity = major / minor diameter, medians. Women's LVOT is more oval at every age (difference constant with age). The short axis runs from the aorto-mitral side to the septum.");

// ---------------- 9 Echo underestimation ----------------
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Results" });
s.addText("A circle from the PLAX diameter misses ~30% of LVOT area", { placeholder: "title" });
const lev = ["Aortic annulus", "LVOT 2 mm (closed)", "LVOT 4 mm (chord-closed)"];
const levLab = ["Aortic annulus", "2 mm from annulus", "4 mm from annulus"];
const echo = (sx) => lev.map((l) => area.find((r) => r.level === l && r.sex === sx).echo_circle_pct_mean);
s.addChart(pres.charts.BAR, [{ name: "Women", labels: levLab, values: echo("F") }, { name: "Men", labels: levLab, values: echo("M") }],
  { x: 0.6, y: 1.4, w: 7.6, h: 5.4, barDir: "col", barGrouping: "clustered", ...chartFrame(), title: "Echo circle π/4·D² vs measured area",
    chartColors: [HEX.women, HEX.men], showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: '0"%"',
    dataLabelColor: HEX.navy, dataLabelFontSize: 13, dataLabelFontBold: true, valAxisMinVal: -40, valAxisMaxVal: 0,
    valAxisMajorUnit: 10, valAxisLabelFormatCode: '0"%"', barGapWidthPct: 60, objectName: "Echo underestimation chart" });
card(s, 8.6, 1.5, 4.1, 5.2, "Echo explanation card");
s.addText([
  { text: "Why", options: { bold: true, fontSize: 18, color: C.text2, breakLine: true } },
  { text: "The PLAX view cuts the LVOT about 14° from its short axis", options: { bullet: true, breakLine: true } },
  { text: "Squaring that diameter ignores the long side-to-side axis", options: { bullet: true, breakLine: true } },
  { text: " ", options: { fontSize: 8, breakLine: true } },
  { text: "So what", options: { bold: true, fontSize: 18, color: C.text2, breakLine: true } },
  { text: "Stroke volume and continuity-equation AVA scale with this area", options: { bullet: true, breakLine: true } },
  { text: "Women are underestimated by ~4–6 points more than men", options: { bullet: true, breakLine: true } },
  { text: "The annulus, being rounder, is underestimated far less (−9%)", options: { bullet: true } },
], { x: 8.85, y: 1.7, w: 3.65, h: 4.8, fontSize: 14, color: C.text1, valign: "top", paraSpaceAfter: 6, isTextBox: true, objectName: "Echo explanation text" });
s.addNotes("Echo-style PLAX diameter = edge-to-edge line through the outline's centre, in the plane containing the aortic axis and the mitral centre. End-diastolic geometry; clinical LVOT diameters are measured in mid-systole.");

// ---------------- 10 Area methods ----------------
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Results" });
s.addText("Any ellipse model gets within ~6% of the true area", { placeholder: "title" });
s.addImage({ path: FIG("analysis/fig_area_methods_demo.png"), x: 0.6, y: 1.3, w: 9.86, h: 5.6, objectName: "Area methods figure" });
s.addText([
  { text: "Ellipse methods", options: { bold: true, color: C.text2, breakLine: true } },
  { text: "Moment +1%", options: { breakLine: true } },
  { text: "Feret −2 to +6%", options: { breakLine: true } },
  { text: "Minor–major −4 to +5%", options: { breakLine: true } },
  { text: " ", options: { fontSize: 8, breakLine: true } },
  { text: "Echo circle", options: { bold: true, color: C.accent1, breakLine: true } },
  { text: "−9% annulus", options: { breakLine: true } },
  { text: "−29 to −36% LVOT", options: {} },
], { x: 10.7, y: 1.5, w: 2.1, h: 4.8, fontSize: 14, color: C.text1, valign: "top", isTextBox: true, objectName: "Area methods summary" });
s.addNotes("Top: one LVOT outline (2 mm) with each method's diameters and assumed shape. Bottom: mean ± SD % difference from the measured outline area across all meshes, women lighter, men darker.");

// ---------------- 11 Caveats ----------------
pres.addSection({ title: "Wrap-up" });
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Wrap-up" });
s.addText("What these numbers can and cannot say", { placeholder: "title" });
const CAV = [
  ["End-diastole only", "No systolic meshes exist in this dataset; clinical LVOT diameters are mid-systolic"],
  ["No valve leaflets", "The anterior mitral leaflet is missing, so the closed tube is short and 4 mm uses a chord"],
  ["Averaged hearts", "Each mesh averages many people: results describe typical hearts, not individual spread"],
  ["Annulus = cap contact", "Where the aortic cap meets the LV wall; may sit inside the true hinge points"],
  ["Cross-sectional", "Age differences are between groups, not changes within a person"],
];
CAV.forEach(([h, d], i) => {
  const y = 1.45 + i * 1.06;
  card(s, 0.6, y, 12.1, 0.9, `Caveat ${i + 1} card`);
  numberDot(s, 0.85, y + 0.2, i + 1, C.accent1, `Caveat ${i + 1} dot`);
  s.addText(h, { x: 1.55, y: y + 0.1, w: 3.2, h: 0.7, fontSize: 17, bold: true, color: C.text2, margin: 0, valign: "middle", isTextBox: true, objectName: `Caveat ${i + 1} header` });
  s.addText(d, { x: 4.8, y: y + 0.1, w: 7.7, h: 0.7, fontSize: 15, color: C.text1, margin: 0, valign: "middle", isTextBox: true, objectName: `Caveat ${i + 1} detail` });
});
s.addNotes("Comparisons by sex and age are more reliable than absolute values. Confidence intervals and p-values are descriptive because participant counts per mesh are unknown.");

// ---------------- 12 Reproducibility and next steps ----------------
s = pres.addSlide({ masterName: "Title only", sectionTitle: "Wrap-up" });
s.addText("Fully reproducible, with clear next steps", { placeholder: "title" });
card(s, 0.6, 1.5, 5.9, 5.2, "Reproducibility card");
s.addText([
  { text: "Reproducible pipeline", options: { bold: true, fontSize: 20, color: C.text2, breakLine: true } },
  { text: "One command rebuilds everything (run_all.sh)", options: { bullet: true, breakLine: true } },
  { text: "All parameters in config.yaml; pinned package versions", options: { bullet: true, breakLine: true } },
  { text: "Every step documented with its code (PIPELINE.md, flowchart)", options: { bullet: true, breakLine: true } },
  { text: "Decisions log with reasons and evidence", options: { bullet: true, breakLine: true } },
  { text: "Run log with code and config fingerprints; git history", options: { bullet: true } },
], { x: 0.9, y: 1.75, w: 5.3, h: 4.8, fontSize: 15, color: C.text1, valign: "top", paraSpaceAfter: 8, isTextBox: true, objectName: "Reproducibility text" });
card(s, 6.8, 1.5, 5.9, 5.2, "Next steps card");
s.addText([
  { text: "Next steps", options: { bold: true, fontSize: 20, color: C.accent1, breakLine: true } },
  { text: "Calibrate an echo LVOT-area correction on paired systolic echo and CT/3D data", options: { bullet: true, breakLine: true } },
  { text: "Test corrected AVA against outcomes; thresholds may need to change", options: { bullet: true, breakLine: true } },
  { text: "Make any correction sex- and level-specific", options: { bullet: true, breakLine: true } },
  { text: "Revisit PLAX-plane sensitivity when the echo comparison is developed", options: { bullet: true } },
], { x: 7.1, y: 1.75, w: 5.3, h: 4.8, fontSize: 15, color: C.text1, valign: "top", paraSpaceAfter: 8, isTextBox: true, objectName: "Next steps text" });
s.addNotes("Code and docs: ~/Documents/Claude/lvot-ukb (README, WRITEUP, docs/PIPELINE.md, docs/DECISIONS.md).");

const out = path.join(ROOT, "outputs/LVOT_highlights.pptx");
pres.writeFile({ fileName: out }).then(async () => {
  const applyThemePath = process.env.APPLY_THEME;
  if (applyThemePath) {
    const { applyTheme } = require(applyThemePath);
    await applyTheme(out, THEME);
  }
  console.log("wrote", out);
});
