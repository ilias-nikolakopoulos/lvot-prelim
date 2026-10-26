// STEP 7d — One-slide deck: "What is a mesh?" (full-slide figure from scripts/16_mesh_explainer.py + speaker notes).
// Usage: NODE_PATH=<dir with pptxgenjs> node scripts/17_mesh_slide.js   -> outputs/slide_what_is_a_mesh.pptx
const path = require("path");
const pptxgen = require("pptxgenjs");
const ROOT = path.resolve(__dirname, "..");
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.333 x 7.5 in, same 16:9 as the figure
pres.title = "What is a mesh?";
pres.author = "Ilias Nikolakopoulos";
const s = pres.addSlide();
s.background = { color: "FFFFFF" };
s.addImage({ path: path.join(ROOT, "outputs/figures/mesh_explainer.png"), x: 0, y: 0, w: 13.333, h: 7.5,
  altText: "What is a mesh: the whole heart coloured by region, a zoom showing surface triangles, a slice showing the wall filled with tetrahedra, and one tetrahedron with its four nodes",
  objectName: "Mesh explainer figure" });
s.addNotes([
  "A mesh is a digital 3D model of the heart built from about 2.1 million small tetrahedra (four-cornered solids), joined at about 430,000 shared corner points.",
  "Each tetrahedron carries a region tag: LV wall, RV wall, or one of four flat caps that close the mitral, aortic, tricuspid and pulmonary openings. The caps are not real valves; there are no leaflets.",
  "Each corner point has x, y, z coordinates (stored in micrometres, converted to mm).",
  "These meshes are averages of UK Biobank participants with the same sex, age and BMI, at end-diastole.",
  "Every measurement in this project is made directly on these points and triangles.",
].join("\n"));
pres.writeFile({ fileName: path.join(ROOT, "outputs/slide_what_is_a_mesh.pptx") }).then((f) => console.log("wrote", f));
