// STEP 7f — Ellipticity slide (dark green, Roboto), two messages.
//   outputs/slide_ellipticity_static.pptx    two slides: build 1, then build 2
//   outputs/slide_ellipticity_animated.pptx  one slide; the on-click build is added by scripts/20_add_click_build.py
// Images: outputs/figures/slides/ellipticity_depth_dark_build{1,2}.png (scripts/19_ellipticity_depth_slide.py).
// Usage: NODE_PATH=<dir with pptxgenjs> node scripts/20_ellipticity_slide.js
const path = require("path");
const pptxgen = require("pptxgenjs");
const ROOT = path.resolve(__dirname, "..");
const IMG = (b) => path.join(ROOT, `outputs/figures/slides/ellipticity_depth_dark_build${b}.png`);

const BG = "00301F", INK = "F5F0E6", WOMEN = "D55E00", MUTED = "9FB3A8", FONT = "Roboto";
const L = 0.6;                         // shared left edge for titles, chart axis and footnote
const CH = { x: 0.42, y: 1.75, w: 11.2, h: 11.2 * 1240 / 2666 };  // chart keeps the PNG aspect ratio

function line1(s) {
  s.addText("The LVOT becomes more oval below the annulus", {
    x: L, y: 0.38, w: 12.1, h: 0.72, fontFace: FONT, fontSize: 38, bold: true, color: INK, margin: 0,
    valign: "middle", isTextBox: true, objectName: "Headline depth" });
}
function line2(s) {
  s.addText([
    { text: "…and ", options: { color: INK } },
    { text: "women's", options: { color: WOMEN, bold: true } },
    { text: " LVOT is more oval than men's", options: { color: INK } },
  ], { x: L, y: 1.08, w: 12.1, h: 0.58, fontFace: FONT, fontSize: 30, margin: 0, valign: "middle",
       isTextBox: true, objectName: "Headline sex" });
}
function chart(s, b) {
  s.addImage({ path: IMG(b), ...CH, objectName: `Chart build ${b}`,
    altText: b === 1
      ? "Box plots of LVOT ellipticity at the annulus, 2 mm and 4 mm below it, rising from 1.14 to about 1.5"
      : "Same box plots coloured by sex: no sex gap at the annulus, women more oval by 0.07 at 2 mm and 0.10 at 4 mm" });
}
function footnote(s) {
  s.addText("Ellipticity = major / minor diameter (minor–major axis ellipse) · 622 matched same-age, same-BMI pairs · UK Biobank meshes, end-diastole",
    { x: L, y: 7.05, w: 12.1, h: 0.3, fontFace: FONT, fontSize: 12, color: MUTED, margin: 0, isTextBox: true,
      objectName: "Footnote" });
}
const NOTES1 = "The LVOT is nearly round at the annulus (1.14) and becomes clearly oval within a few millimetres below it: about 1.4–1.5 at 2 mm and 1.43–1.53 at 4 mm.";
const NOTES2 = "No sex difference at the annulus; women's LVOT is more oval than men's by 0.07 at 2 mm and 0.10 at 4 mm, at every age. Matched same-age, same-BMI pairs.";

function deck() {
  const p = new pptxgen();
  p.layout = "LAYOUT_WIDE";
  p.author = "Ilias Nikolakopoulos";
  p.title = "LVOT ellipticity by depth and sex";
  return p;
}

// static: two slides
const st = deck();
for (const b of [1, 2]) {
  const s = st.addSlide();
  s.background = { color: BG };
  line1(s);
  if (b === 2) line2(s);
  chart(s, b);
  footnote(s);
  s.addNotes(b === 1 ? NOTES1 : NOTES2);
}
// animated: one slide with everything; build-2 objects are hidden until the click (added afterwards)
const an = deck();
const s = an.addSlide();
s.background = { color: BG };
line1(s);
line2(s);
chart(s, 1);
chart(s, 2);
footnote(s);
s.addNotes(NOTES1 + "\n[click] " + NOTES2);

Promise.all([
  st.writeFile({ fileName: path.join(ROOT, "outputs/slide_ellipticity_static.pptx") }),
  an.writeFile({ fileName: path.join(ROOT, "outputs/slide_ellipticity_animated.pptx") }),
]).then((f) => console.log("wrote", f.join(", ")));
