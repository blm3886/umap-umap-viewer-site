import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import Papa from 'papaparse';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';

const FILE_PATH = './saturn_joint_umap_3d_1.csv';
//spawnFloatingGeneNamesForSet(setName);
//clearFloatingGeneNames();

// === PAGE STYLE ===
const __style = document.createElement('style');
__style.textContent = `
  html, body { margin: 0; padding: 0; }
  canvas { display: block; }
`;
document.head.appendChild(__style);

// === FLOATING GENE NAME STYLES (high-contrast) ===
const __floatStyle = document.createElement('style');
__floatStyle.textContent = `
  .gene-float {
    position: absolute;
    pointer-events: none;
    z-index: 2000; /* sit above canvas */
    color: #ffffff;
    font: 800 30px/1.3 system-ui, -apple-system, Segoe UI, Arial;
    white-space: nowrap;

    /* Strong outline for any background */
    -webkit-text-stroke: 1px rgba(0,0,0,0.55);
    text-shadow:
      0 1px 0 rgba(0,0,0,.35),
      0 0 4px rgba(0,0,0,.55),
      0 0 10px rgba(0,0,0,.65);

    /* Soft drift */
    animation: geneFloatDrift linear infinite;
  }
  .gene-float::before {
    /* Subtle pill backdrop to boost contrast without feeling like a panel */
    content: "";
    position: absolute;
    inset: -3px -8px;           /* padding around text */
    background: rgba(0,0,0,0.48);
    border-radius: 10px;
    filter: blur(0.5px);
    z-index: -1;
  }
  .gene-float.magic {
    color: #00e5ff;             /* neon cyan for "magic" genes */
    -webkit-text-stroke: 1px rgba(0,20,30,0.65);
    text-shadow:
      0 0 6px rgba(0,229,255,0.65),
      0 0 16px rgba(0,229,255,0.45),
      0 1px 0 rgba(0,0,0,.45);
  }

  @keyframes geneFloatDrift {
    0%   { transform: translate(0px, 0px);   opacity:.98; }
    50%  { transform: translate(10px, -8px); opacity:.92; }
    100% { transform: translate(0px, 0px);   opacity:.98; }
  }
`;
document.head.appendChild(__floatStyle);

const MAGIC_GENES = new Set([
  "Tarv.1014.HPI3.6.g238150",
  "LotjaGi5g1v0093600",
  "Tarv.1014.HPI3.1.g047610",
  "Gmax.la3.HPI3.Chr02.g031820",
  "Tarv.1014.HPI3.1.g047580",
  "Tarv.1014.HPI3.3.g114450",
  "AT4G01130",
]);

function normalizeId(id) {
  if (id == null) return '';
  let s = String(id).trim();
  // IMPORTANT: do NOT strip "-1-0" if your expr files include it.
  // If you had `s = s.replace(/-\d+-\d+$/, '');` remove/comment it.
  s = s.replace(/^"+|"+$/g, ''); // only strip accidental quotes
  return s;
}


// === LAYOUT ===
let LEFT_PANE_PX = Math.round(window.innerWidth * 0.22);
const PANE_GAP_PX = 8;
const PANE_INNER_PAD = 6;

// const AI_PRED_SETS = {
//   MACROGENE_SET1: [
//     "AT1G01010","AT2G22222","Gmax.14G077770","Tarv.00001","AT4G12345"
//   ],
//   MACROGENE_SET2: [
//     "AT3G33333","AT5G88888","AT1G05050","Gmax.01G123400"
//   ]
// };


// === SPECIES COLORS ===
const speciesColors = {
  arabidopsis: "#00FFD1",
  pennycress:  "#FF61F6",
  soybean:     "#7CFC00",
  rice:        "#F4A261",
  lotus:       "#FFD700",
  sorghum:     "#5A5DFF"
};

const BLOOM_NORMAL = {
  strength: 1.0,  
  radius:   0.3,
  threshold:0.6,
  exposure: 1.6,
  hlSize:   0.085, 
  hlOpacity:0.95
};

const BLOOM_FOR_SET = {
  strength: 0.22, 
  radius:   0.25,
  threshold:0.85,
  exposure: 1.2,
  hlSize:   0.05,
  hlOpacity:0.35
};

function applyGlowPreset(p) {
  bloomPass.strength  = p.strength;
  bloomPass.radius    = p.radius;
  bloomPass.threshold = p.threshold;

  renderer.toneMappingExposure = p.exposure;

  if (highlightPoints && highlightPoints.material) {
    highlightPoints.material.size    = p.hlSize;
    highlightPoints.material.opacity = p.hlOpacity;
    highlightPoints.material.needsUpdate = true;
  }
}

// GENES = [
//     # "",
//     # "LotjaGi5g1v0093600",
//     # "Tarv.1014.HPI3.1.g047610",
//     "Gmax.la3.HPI3.Chr02.g031820",
//     # "Tarv.1014.HPI3.1.g047580",
//     # "Tarv.1014.HPI3.3.g114450",
//     # "AT1G28600",
//     # "AT4G01130",
//     # "AT1G28610"
// ]

// === GENE SETS ===
// Put the exact gene identifiers you have CSVs for (expr_<ID>.csv).
// Example: expr_AT1G01010.csv, expr_Gmax.01G123400.csv, etc.
const GENE_SETS = {
  MACROGENE_SET_1: ["Tarv.1014.HPI3.6.g238150", "LotjaGi5g1v0093600","Tarv.1014.HPI3.1.g047610","Gmax.la3.HPI3.Chr02.g031820","Tarv.1014.HPI3.1.g047580","Tarv.1014.HPI3.3.g114450","AT1G28600_GGL2_GDSL(HPI)","AT4G01130", "AT1G28610_GGL3_GDSL(HPI)"],
  MACROGENE_SET_2: ["AT3G22222", "Gmax.14G077770"]
};
// "max" lights a cell if any gene in the set is expressed; "sum" stacks contributions.
const DEFAULT_SET_AGG = "max";
let AUTO_SHOW_ONLY_ON_SET = true;

// === STATE: main + overlay ===
let cellIndex = new Map();     
let exprMain = null;           
let highlightPoints = null;    
let currentMin = 0, currentMax = 1;
let baseColors = null;        
let activeGene = null;      

// === THUMBNAIL RIGHT PANE (one per species) ===
const PLANT_PANELS = [
  { label: 'Arabidopsis', filter: (r) => (r.species || '').toLowerCase() === 'arabidopsis' },
  { label: 'Pennycress',  filter: (r) => (r.species || '').toLowerCase() === 'pennycress' },
  { label: 'Soybean',     filter: (r) => (r.species || '').toLowerCase() === 'soybean' },
  { label: 'Rice',        filter: (r) => (r.species || '').toLowerCase() === 'rice' },
  { label: 'Lotus',       filter: (r) => (r.species || '').toLowerCase() === 'lotus' },
  { label: 'Sorghum',     filter: (r) => (r.species || '').toLowerCase() === 'sorghum' },
];
const NUM_THUMBS = PLANT_PANELS.length;

// === SCENE: MAIN (cross-species UMAP) ===
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x000000);

// Camera: MAIN
const camera = new THREE.PerspectiveCamera(75, 1, 0.1, 2000);
camera.position.set(0, 0, 20);

// Renderer
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setScissorTest(true);
renderer.outputColorSpace = THREE.SRGBColorSpace;
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.toneMappingExposure = 1.6;
document.body.appendChild(renderer.domElement);

// Post FX
const composer = new EffectComposer(renderer);
const renderScene = new RenderPass(scene, camera);
const bloomPass = new UnrealBloomPass(
  new THREE.Vector2(window.innerWidth, window.innerHeight),
  1.0,  
  0.3,  
  0.6  
);
composer.addPass(renderScene);
composer.addPass(bloomPass);

// Controls (MAIN)
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;

// Lights (MAIN)
scene.add(new THREE.AmbientLight(0xffffff, 0.4));
const keyLight = new THREE.PointLight(0xffffff, 1.4);
keyLight.position.set(20, 20, 20);
scene.add(keyLight);

// === THUMBNAIL SCENES ===
const thumbs = []; // { scene, camera, points, label }
for (let i = 0; i < NUM_THUMBS; i++) {
  const scn = new THREE.Scene();
  scn.background = new THREE.Color(0x0b0b0e);
  const cam = new THREE.PerspectiveCamera(55, 1, 0.1, 2000);
  cam.position.set(0, 0, 25);
  scn.add(new THREE.AmbientLight(0xffffff, 0.5));
  scn.add(new THREE.DirectionalLight(0xffffff, 0.8));
  thumbs.push({ scene: scn, camera: cam, points: null, label: PLANT_PANELS[i].label });
}

// === DATA LOADING FOR MAIN SCREEN ===
let centeredData = null;
Papa.parse(FILE_PATH, {
  download: true,
  header: true,
  dynamicTyping: true,
  complete: function(result) {
    centeredData = preprocessData(result.data);

    // build cell_id -> point index
    centeredData.forEach((r, idx) => {
    const nid = normalizeId(r.cell_id ?? r.plant_id ?? `row_${idx}`);
    if (nid) cellIndex.set(nid, idx);
    });

    //centeredData.forEach((r, idx) => { if (r.cell_id) cellIndex.set(String(r.cell_id), idx); });

    plotMain(centeredData);
    plotThumbnails(centeredData);
    applyMainCameraAspect();
    layoutLabels();
    injectGeneUI();
  },
  error: function(err) {
    console.error("Error loading CSV:", err);
  }
});

// === CENTER & PREPROCESS ===
function preprocessData(data) {
  const rows = data.filter(r => Number.isFinite(r.x) && Number.isFinite(r.y) && Number.isFinite(r.z));

  const mean = rows.reduce((acc, r) => { acc.x += r.x; acc.y += r.y; acc.z += r.z; return acc; }, {x:0,y:0,z:0});
  mean.x /= rows.length; mean.y /= rows.length; mean.z /= rows.length;

  return rows.map((r, i) => {
    const species = (r.species || '').toLowerCase();
    const colorHex = (r.ai_flag === 1 || r.ai_flag === "1") ? "#ff33cc" : (speciesColors[species] || "cyan");
    const c = new THREE.Color(colorHex);
    const cid = (r.cell_id ?? r.plant_id ?? `row_${i}`) + '';

    return {
      x: r.x - mean.x,
      y: r.y - mean.y,
      z: r.z - mean.z,
      species: r.species,
      label: r.label,
      cell_id: cid,
      color: [c.r, c.g, c.b]
    };
  });
}

// === MAIN POINTS ===
let mainPoints = null;
function plotMain(prepped) {
  const N = prepped.length;
  const positions = new Float32Array(N * 3);
  const colors    = new Float32Array(N * 3);
  exprMain        = new Float32Array(N); 

  for (let i = 0, k = 0; i < N; i++, k += 3) {
    const r = prepped[i];
    positions[k]   = r.x; positions[k+1] = r.y; positions[k+2] = r.z;
    colors[k]      = r.color[0]; colors[k+1] = r.color[1]; colors[k+2] = r.color[2];
  }
  baseColors = new Float32Array(colors);

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("color",    new THREE.BufferAttribute(colors,    3));
  geometry.setAttribute("expr",     new THREE.BufferAttribute(exprMain,  1));

  const mat = new THREE.PointsMaterial({
    size: 0.05,
    sizeAttenuation: true,
    vertexColors: true,
    transparent: true,
    opacity: 0.9,
    depthWrite: false,
    blending: THREE.AdditiveBlending
  });

  if (mainPoints) scene.remove(mainPoints);
  mainPoints = new THREE.Points(geometry, mat);
  scene.add(mainPoints);

  // overlay for picked/highlighted cells
  if (!highlightPoints) {
    const hg = new THREE.BufferGeometry();
    hg.setAttribute("position", new THREE.Float32BufferAttribute([], 3));
    const hmat = new THREE.PointsMaterial({
      size: 0.085,
      sizeAttenuation: true,
      color: 0xffff88,
      transparent: true,
      opacity: 0.95,
      blending: THREE.AdditiveBlending,
      depthWrite: false
    });
    highlightPoints = new THREE.Points(hg, hmat);
    scene.add(highlightPoints);
  }
}

// === THUMB POINTS ===
function plotThumbnails(prepped) {
  for (let i = 0; i < NUM_THUMBS; i++) {
    const cfg = PLANT_PANELS[i];
    const subset = prepped.filter(cfg.filter);
    const points = buildPoints(subset);

    if (thumbs[i].points) thumbs[i].scene.remove(thumbs[i].points);
    thumbs[i].points = points;
    thumbs[i].scene.add(points);

    fitCameraToPoints(thumbs[i].camera, points, 1.4);
    ensureLabel(i, thumbs[i].label);
  }
}

function buildPoints(subset) {
  const n = subset.length;
  const positions = new Float32Array(n * 3);
  const colors = new Float32Array(n * 3);

  for (let i = 0; i < n; i++) {
    const r = subset[i];
    const ix = i * 3;
    positions[ix]     = r.x;
    positions[ix + 1] = r.y;
    positions[ix + 2] = r.z;
    colors[ix]     = r.color[0];
    colors[ix + 1] = r.color[1];
    colors[ix + 2] = r.color[2];
  }

  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  geom.setAttribute('color', new THREE.BufferAttribute(colors, 3));

  return new THREE.Points(geom, new THREE.PointsMaterial({
    size: 0.045,
    sizeAttenuation: true,
    vertexColors: true,
    transparent: false,
    opacity: 1.0,
    depthWrite: true,
    blending: THREE.NormalBlending
  }));
}

// === CAMERA HELPERS ===
function fitCameraToPoints(camera, points, pad = 1.3) {
  const box = new THREE.Box3().setFromObject(points);
  const size = new THREE.Vector3(); box.getSize(size);
  const center = new THREE.Vector3(); box.getCenter(center);
  const maxDim = Math.max(size.x, size.y, size.z);
  const fov = camera.fov * (Math.PI / 180);
  let dist = (maxDim / 2) / Math.tan(fov / 2);
  dist *= pad;

  camera.position.set(center.x, center.y, dist + center.z);
  camera.lookAt(center);
  camera.updateProjectionMatrix();
}

// === LABELS (for thumbnails) ===
function ensureLabel(i, text) {
  const id = `thumb-label-${i}`;
  let el = document.getElementById(id);
  if (!el) {
    el = document.createElement('div');
    el.id = id;
    el.style.position = 'absolute';
    el.style.pointerEvents = 'none';
    el.style.textAlign = 'center';
    el.style.font = '700 15px Arial, sans-serif';
    el.style.color = 'silver';
    el.style.textShadow = '0 0 2px #999, 0 0 6px #ccc';
    document.body.appendChild(el);
  }
  el.textContent = text;
}

function layoutLabels() {
  for (let i = 0; i < NUM_THUMBS; i++) {
    const r = getThumbRect(i);
    const el = document.getElementById(`thumb-label-${i}`);
    if (!el) continue;
    el.style.left = r.x + 'px';
    el.style.top  = (r.y + 6) + 'px';
    el.style.width = r.w + 'px';
  }
}

// === LAYOUT / VIEWPORTS ===
function getThumbRect(index) {
  const totalHeight = window.innerHeight;
  const eachH = Math.floor((totalHeight - (NUM_THUMBS + 1) * PANE_GAP_PX) / NUM_THUMBS);
  const x = 0;
  const y = PANE_GAP_PX + index * (eachH + PANE_GAP_PX);
  const w = LEFT_PANE_PX;
  const h = eachH;
  return { x, y, w, h };
}

function getMainRect() {
  return {
    x: LEFT_PANE_PX + PANE_GAP_PX,
    y: 0,
    w: window.innerWidth - LEFT_PANE_PX - PANE_GAP_PX,
    h: window.innerHeight
  };
}

function applyMainCameraAspect() {
  const mr = getMainRect();
  camera.aspect = Math.max(0.1, mr.w / mr.h);
  camera.updateProjectionMatrix();
}

// === RESIZE ===
window.addEventListener('resize', () => {
  LEFT_PANE_PX = Math.round(window.innerWidth * 0.22);
  renderer.setSize(window.innerWidth, window.innerHeight);
  composer.setSize(window.innerWidth, window.innerHeight);

  applyMainCameraAspect();

  for (const t of thumbs) {
    t.camera.aspect = 1;
    t.camera.updateProjectionMatrix();
  }
  layoutLabels();
});

// === SINGLE-GENE CSV LOADER (kept for convenience) ===
async function loadExpressionCSV(path){
  const text = await fetch(path).then(r => { if(!r.ok) throw new Error(`failed to load ${path}`); return r.text(); });
  const lines = text.trim().split(/\r?\n/);
  const hasHeader = /cell/i.test(lines[0]) || /expr/i.test(lines[0]);
  const start = hasHeader ? 1 : 0;

  exprMain.fill(0);
  let min = Infinity, max = -Infinity;

  for (let i = start; i < lines.length; i++) {
    if (!lines[i]) continue;
    const [id, valStr] = lines[i].split(/,|\t/);
    const idx = cellIndex.get(String(id));
    if (idx === undefined) continue;
    const v = parseFloat(valStr);
    if (!isFinite(v)) continue;
    exprMain[idx] = v;
    if (v < min) min = v; if (v > max) max = v;
  }
  mainPoints.geometry.getAttribute("expr").needsUpdate = true;

  if (!isFinite(min) || !isFinite(max)) { min = 0; max = 1; }
  currentMin = min; currentMax = max;

  console.log(currentMin)
  console.log(currentMax)
  console.log("__________________________________________________!!!!")

  const thr = currentMin + 0.5 * (currentMax - currentMin);
  recolorByExpression(thr);
  updateHighlight(thr);
  updateThrUI(thr);
}

// === SET HELPERS ===
// function parseExprCSVToMap(text) {
//   const lines = text.trim().split(/\r?\n/);
//   const hasHeader = /cell/i.test(lines[0]) || /expr/i.test(lines[0]);
//   const start = hasHeader ? 1 : 0;

//   const m = new Map();
//   for (let i = start; i < lines.length; i++) {
//     if (!lines[i]) continue;
//     const [id, valStr] = lines[i].split(/,|\t/);
//     const v = parseFloat(valStr);
//     if (!isFinite(v)) continue;
//     m.set(String(id), v);
//   }
//   return m;
// }
function parseExprCSVToMap(text) {
  const lines = text.trim().split(/\r?\n/);
  const hasHeader = /cell/i.test(lines[0]) || /expr/i.test(lines[0]);
  const start = hasHeader ? 1 : 0;
  const m = new Map();
  for (let i = start; i < lines.length; i++) {
    if (!lines[i]) continue;
    const [id, valStr] = lines[i].split(/,|\t/);
    const nid = normalizeId(id);
    const v = parseFloat(valStr);
    if (!nid || !isFinite(v)) continue;
    m.set(nid, v);
  }
  return m;
}


function aggregateExprIntoArray(exprMaps, mode = "max") {
  exprMain.fill(0);
  let min = +Infinity, max = -Infinity;

  for (const [id, idx] of cellIndex.entries()) {
    let val = 0;
    if (mode === "sum") {
      let acc = 0;
      for (const m of exprMaps) acc += (m.get(id) || 0);
      val = acc;
    } else {
      // default to max
      let best = 0;
      for (const m of exprMaps) {
        const v = m.get(id) || 0;
        if (v > best) best = v;
      }
      val = best;
    }
    exprMain[idx] = val;
    if (val < min) min = val;
    if (val > max) max = val;
  }

  mainPoints.geometry.getAttribute("expr").needsUpdate = true;

  if (!isFinite(min) || !isFinite(max)) { min = 0; max = 1; }
  currentMin = min;
  currentMax = max;
}

async function loadExpressionSET(setName, mode = DEFAULT_SET_AGG) {

  const genes = GENE_SETS[setName] || [];
  if (genes.length === 0) {
    console.warn(`GENE_SETS["${setName}"] is empty.`);
    exprMain.fill(0);
    mainPoints.geometry.getAttribute("expr").needsUpdate = true;
    recolorByExpression(0);
    updateThrUI(0);
    return;
  }

  const texts = await Promise.all(
    genes.map(g =>
      fetch(`expr_${g}.csv`).then(r => {
        if (!r.ok) throw new Error(`failed to load expr_${g}.csv`);
        return r.text();
      })
    )
  );

  const maps = texts.map(parseExprCSVToMap);
  aggregateExprIntoArray(maps, mode);

  const thr = currentMin + 0.1 * (currentMax - currentMin);
  //const thr = currentMin
  recolorByExpression(thr);
  updateHighlight(thr);
  updateThrUI(thr);
}

// === MODES: show-only or recolor ===
function showOnlyExpressing(threshold){
  const geom  = mainPoints.geometry;
  const posA  = geom.getAttribute("position").array;
  const exprA = geom.getAttribute("expr").array;

  const picked = [];
  for (let i = 0; i < exprA.length; i++) {
    if (exprA[i] >= threshold) {
      const k = i * 3;
      picked.push(posA[k], posA[k+1], posA[k+2]);
    }
  }
  const hg = highlightPoints.geometry;
  hg.setAttribute("position", new THREE.Float32BufferAttribute(picked, 3));
  hg.attributes.position.needsUpdate = true;
  hg.computeBoundingSphere();

  mainPoints.visible = false;
  highlightPoints.visible = picked.length > 0;
}

function showAllCells(){
  mainPoints.visible = true;
  highlightPoints.visible = false;
}

function resetColorsToSpecies(){
  if (!baseColors) return;
  const colors = mainPoints.geometry.getAttribute("color").array;
  colors.set(baseColors);
  mainPoints.geometry.getAttribute("color").needsUpdate = true;
}

function recolorByExpression(threshold){
  const geom     = mainPoints.geometry;
  const colors   = geom.getAttribute("color").array;
  const exprA    = geom.getAttribute("expr").array;
  if (!exprA || !baseColors) return;

  const GRAY = 0.5;
  for (let i = 0; i < exprA.length; i++) {
    const k = i * 3;
    if (exprA[i] >= threshold) {
      colors[k]   = baseColors[k];
      colors[k+1] = baseColors[k+1];
      colors[k+2] = baseColors[k+2];
    } else {
      colors[k]   = GRAY;
      colors[k+1] = GRAY;
      colors[k+2] = GRAY;
    }
  }
  geom.getAttribute("color").needsUpdate = true;
}

function updateHighlight(threshold){
  if (mainPoints) highlightPoints.rotation.copy(mainPoints.rotation);
  const geom  = mainPoints.geometry;
  const posA  = geom.getAttribute("position").array;
  const exprA = geom.getAttribute("expr").array;

  const picked = [];
  for (let i = 0; i < exprA.length; i++) {
    if (exprA[i] >= threshold) {
      const k = i*3;
      picked.push(posA[k], posA[k+1], posA[k+2]);
    }
  }
  const hg = highlightPoints.geometry;
  hg.setAttribute("position", new THREE.Float32BufferAttribute(picked, 3));
  hg.attributes.position.needsUpdate = true;
  hg.computeBoundingSphere();
  highlightPoints.visible = picked.length > 0;
}

// === FLOATING GENE NAMES (sprinkled over main viewport) ===
let floatingGeneEls = [];

function clearFloatingGeneNames() {
  for (const el of floatingGeneEls) if (el?.parentNode) el.parentNode.removeChild(el);
  floatingGeneEls = [];
}

// Random position *outside* a central “safe” box to keep labels readable
function randomPositionOutsideCenter(mr, insetPx) {
  // inner box we avoid
  const inner = {
    x: mr.x + insetPx,
    y: mr.y + insetPx,
    w: Math.max(1, mr.w - insetPx * 2),
    h: Math.max(1, mr.h - insetPx * 2)
  };

  // rejection sample: try up to N times to avoid the inner box
  for (let tries = 0; tries < 30; tries++) {
    const x = mr.x + Math.random() * mr.w;
    const y = mr.y + Math.random() * mr.h;

    const insideInner =
      x >= inner.x && x <= inner.x + inner.w &&
      y >= inner.y && y <= inner.y + inner.h;

    if (!insideInner) return { x, y };
  }

  return { x: mr.x + Math.random() * mr.w, y: (Math.random() < 0.5 ? mr.y + 6 : mr.y + mr.h - 6) };
}

function spawnFloatingGeneNamesForSet(setName) {
  //clearFloatingGeneNames();
  const genes = GENE_SETS[setName] || [];
  if (!genes.length) return;

  const mr = getMainRect();
  const inset = Math.max(80, Math.min(mr.w, mr.h) * 0.25);

  genes.forEach((g) => {
  let { x, y } = randomPositionOutsideCenter(mr, inset);

  const el = document.createElement('span');
  el.className = 'gene-float' + (MAGIC_GENES.has(g) ? ' magic' : '');
  el.textContent = MAGIC_GENES.has(g) ? `✨ ${g}` : g;

  document.body.appendChild(el);
  floatingGeneEls.push(el);

  // measure element width/height after inserting
  const rect = el.getBoundingClientRect();
  const w = rect.width;
  const h = rect.height;

  // clamp so it stays within main rect
  const minX = mr.x + 4;
  const maxX = mr.x + mr.w - w - 4;
  const minY = mr.y + 4;
  const maxY = mr.y + mr.h - h - 4;

  x = Math.min(maxX, Math.max(minX, x));
  y = Math.min(maxY, Math.max(minY, y));

  el.style.left = `${Math.round(x)}px`;
  el.style.top  = `${Math.round(y)}px`;

  const dur = 6 + Math.random() * 7;
  el.style.animationDuration = `${dur}s`;
  el.style.animationDelay = `${-Math.random() * dur}s`;
});

}

window.addEventListener('resize', () => {
  if (floatingGeneEls.length && activeGene) {
    //spawnFloatingGeneNamesForSet(activeGene);
  }
});


// let floatingGeneEls = [];

// function clearFloatingGeneNames() {
//   for (const el of floatingGeneEls) {
//     if (el && el.parentNode) el.parentNode.removeChild(el);
//   }
//   floatingGeneEls = [];
// }

// function spawnFloatingGeneNamesForSet(setName) {
//   clearFloatingGeneNames();
//   const genes = GENE_SETS[setName] || [];
//   if (!genes.length) return;

//   const mr = getMainRect(); // {x,y,w,h}

//   genes.forEach((g) => {
//     // random position inside the main viewport
//     const x = mr.x + Math.random() * mr.w;
//     const y = mr.y + Math.random() * mr.h;

//     const el = document.createElement('span');
//     el.className = 'gene-float';
//     el.textContent = `✨ ${g}`;
//     el.style.left = `${Math.round(x)}px`;
//     el.style.top  = `${Math.round(y)}px`;

//     // vary duration + phase so they don't move in sync
//     const dur = 6 + Math.random() * 7;         // 6–13s
//     const delay = (-Math.random() * dur) + 's'; // negative = desync immediately
//     el.style.animationDuration = `${dur}s`;
//     el.style.animationDelay = delay;

//     document.body.appendChild(el);
//     floatingGeneEls.push(el);
//   });
// }

// window.addEventListener('resize', () => {
//   if (floatingGeneEls.length && activeGene) {
//     spawnFloatingGeneNamesForSet(activeGene);
//   }
// });


// === RIGHT PANE UI ===
// <div style="display:flex;gap:6px;align-items:center;margin:8px 0">
//   <select id="aggMode" style="flex:1;background:#1b1f2a;color:#ddd;border:1px solid #333;border-radius:6px;padding:4px">
//     <option value="max" selected>max (any gene lights up)</option>
//     <option value="sum">sum (stack contributions)</option>
//   </select>
// </div>
function injectGeneUI(){
  const ui = document.createElement('div');
  //ui.style.cssText = 'position:absolute;right:14px;top:14px;width:260px;max-height:70vh;overflow:auto;padding:10px 12px;background:rgba(20,20,28,.75);backdrop-filter:blur(4px);border:1px solid #333;border-radius:12px;color:#ddd;font:13px system-ui;z-index:1000';
  ui.style.cssText = 'position:absolute;right:14px;top:14px;width:460px;max-height:70vh;overflow:auto;padding:10px 12px;background:rgba(20,20,28,.75);backdrop-filter:blur(4px);border:1px solid #333;border-radius:12px;color:#ddd;font:13px system-ui;z-index:1000';

  //   ui.innerHTML = `
//     <div style="font-weight:700;margin-bottom:8px">PREDICTED GENE CELL EXPRESSION</div>
//     <div id="geneList"></div>
//     <div style="height:10px"></div>

//     <div style="font-weight:700;margin:6px 0 4px">Threshold</div>
//     <div style="display:flex;gap:8px;align-items:center">
//       <input id="thr" type="range" min="0" max="1" step="0.001" value="0.7" style="flex:1">
//       <span id="thrVal" style="width:58px;text-align:right;color:#aaa">auto</span>
//     </div>

//     <label style="display:flex;gap:6px;align-items:center;margin-top:6px;opacity:.95">
//       <input id="showOnlyChk" type="checkbox">
//       Show only expressing cells
//     </label>
//   `;
    ui.innerHTML = `
    <div id="pgceHeader" style="
        display:flex;align-items:center;justify-content:space-between;
        cursor:pointer;user-select:none;padding:6px 4px;margin:-6px -4px 6px;
        border-radius:8px;
        width:460px;
    ">
        <div style="font-weight:700">PREDICTED GENE CELL EXPRESSION</div>
        <span id="pgceCaret" style="
        display:inline-block;transform:rotate(0deg);transition:transform .18s ease;
        font-weight:900;font-size:16px;color:#bbb;
        ">▾</span>
    </div>

    <div id="pgceBody">
        <div id="geneList"></div>
        <div style="height:10px"></div>

        <div style="font-weight:700;margin:6px 0 4px">Threshold</div>
        <div style="display:flex;gap:8px;align-items:center">
        <input id="thr" type="range" min="0" max="1" step="0.001" value="0.7" style="flex:1">
        <span id="thrVal" style="width:58px;text-align:right;color:#aaa">auto</span>
        </div>

        <label style="display:flex;gap:6px;align-items:center;margin-top:6px;opacity:.95">
        <input id="showOnlyChk" type="checkbox">
        Show only expressing cells
        </label>
    </div>
    `;


  document.body.appendChild(ui);

  // --- Dropdown toggle for the PGCE panel ---
const pgceHeader = ui.querySelector('#pgceHeader');
const pgceBody   = ui.querySelector('#pgceBody');
const pgceCaret  = ui.querySelector('#pgceCaret');

let pgceOpen = true; // default open

function reflowBelowPanels() {
  // If you position panels below this one using its height, update them here.
  const stack = document.getElementById('geneStack');
  const ai    = document.getElementById('aiPredPanel');

  // Compute the top edge just below this panel
  const r = ui.getBoundingClientRect();
  const nextTop = Math.round(r.top + r.height + 10);

  if (stack) stack.style.top = `${nextTop}px`;
  if (ai)    ai.style.top    = `${nextTop}px`;
}

function setPgceOpen(open) {
  pgceOpen = open;
  pgceBody.style.display = open ? 'block' : 'none';
  pgceCaret.style.transform = open ? 'rotate(0deg)' : 'rotate(-90deg)';
  reflowBelowPanels();
}

pgceHeader.addEventListener('click', () => setPgceOpen(!pgceOpen));
// call once to normalize layout (in case you later want it closed by default)
setPgceOpen(false);




  const list = ui.querySelector('#geneList');

  // Buttons for each set
  Object.keys(GENE_SETS).forEach((setName) => {
    const b = document.createElement('button');
    b.textContent = setName;
    b.style.cssText =
      'display:block;width:100%;text-align:left;margin:4px 0;' +
      'padding:6px 8px;background:#2b3a55;border:0;border-radius:8px;' +
      'color:#fff;cursor:pointer';

    b.onclick = async () => {
      if (activeGene === setName) {
        // toggle off
        hideGeneStack();
        activeGene = null;
        b.style.background = '#2b3a55';
        resetColorsToSpecies();
        showOnlyChk.checked = false;
        showAllCells();
        updateThrUI(NaN);
        applyGlowPreset(BLOOM_NORMAL);
        //clearFloatingGeneNames();
        return;
      }

      activeGene = setName;
      list.querySelectorAll('button').forEach(btn => { btn.style.background = '#2b3a55'; });
      b.style.background = '#445c7a';

      //const agg = ui.querySelector('#aggMode').value || DEFAULT_SET_AGG;
      //await loadExpressionSET(setName, agg);
      await loadExpressionSET(setName, DEFAULT_SET_AGG);
      showGeneStack(setName);
      if (AUTO_SHOW_ONLY_ON_SET) {
        showOnlyChk.checked = true;
        const thrEl = ui.querySelector('#thr');
        const thr = currentMin + parseFloat(thrEl.value) * (currentMax - currentMin);
        showOnlyExpressing(thr);   // hides base cloud, shows only set points
        // keep glow softened if you applied the preset trick
        }
      applyGlowPreset(BLOOM_FOR_SET);
      //spawnFloatingGeneNamesForSet(setName);
    };

    list.appendChild(b);

    // === AI PREDICTED GENE CANDIDATES panel (under the main panel / under geneStack if visible) ===
const aiPanel = document.createElement('div');
aiPanel.id = 'aiPredPanel';
aiPanel.style.cssText = `
  position:absolute;
  right:14px;
  top:${ui.offsetTop + ui.offsetHeight + 10}px;
  width:460px;
  max-height:40vh;
  overflow:auto;
  padding:10px 12px;
  background:rgba(20,20,28,.75);
  backdrop-filter:blur(4px);
  border:1px solid #333;
  border-radius:12px;
  color:#ddd;
  font:20px system-ui;
  z-index:998;
`;
aiPanel.innerHTML = `
  <div style="font-weight:700;margin-bottom:8px">AI-Driven Gene Predictions (✨)</div>
  <div id="aiBtnList"></div>
  <div style="height:8px"></div>
  <div id="aiGeneList" style="display:grid;gap:4px;font-size:25px;line-height:1.35"></div>
`;
document.body.appendChild(aiPanel);

// keep this panel pinned under the main UI (or under geneStack if it's shown)
function positionAIPredPanel() {
  const rMain = ui.getBoundingClientRect();
  const stackEl = document.getElementById('geneStack');
  const baseTop = (stackEl && stackEl.style.display !== 'none')
    ? (stackEl.getBoundingClientRect().top + stackEl.getBoundingClientRect().height + 10)
    : (rMain.top + rMain.height + 10);
  aiPanel.style.top = `${Math.round(baseTop)}px`;
}
window.addEventListener('resize', positionAIPredPanel);

// render gene list with ✨ on magic genes
function renderAIGeneList(setName) {
  const listEl = document.getElementById('aiGeneList');
  if (!listEl) return;
  const genes = GENE_SETS[setName] || [];
  if (!genes.length) {
    listEl.innerHTML = `<div style="opacity:.7">No predictions for ${setName}.</div>`;
    return;
  }
  listEl.innerHTML = genes.map(g => {
    const isMagic = MAGIC_GENES.has(g);
    const style = isMagic
      ? 'color:#00e5ff;font-weight:700;text-shadow:0 0 6px rgba(0,229,255,.6)'
      : '';
    const badge = isMagic ? '✨ ' : '';
    return `<div style="${style}">${badge}${g}</div>`;
  }).join('');
}

// build the AI set buttons
(() => {
  const btnHost = document.getElementById('aiBtnList');
  if (!btnHost) return;

  Object.keys(GENE_SETS).forEach(setName => {
    const btn = document.createElement('button');
    btn.textContent = setName;
    btn.style.cssText =
      'display:block;width:100%;text-align:left;margin:4px 0;' +
      'padding:6px 8px;background:#2b3a55;border:0;border-radius:8px;' +
      'color:#fff;cursor:pointer';

    btn.onclick = () => {
      // highlight the active AI button
      btnHost.querySelectorAll('button').forEach(b => b.style.background = '#2b3a55');
      btn.style.background = '#445c7a';

      renderAIGeneList(setName);
      positionAIPredPanel();
    };

    btnHost.appendChild(btn);
  });

  // optional: auto-render the first set on load
  const first = Object.keys()[0];
  if (first) {
    renderAIGeneList(first);
    // visually mark first as active
    const firstBtn = btnHost.querySelector('button');
    if (firstBtn) firstBtn.style.background = '#445c7a';
    positionAIPredPanel();
  }
})();


  });

// === floating gene list container (below the panel) ===
    const geneStack = document.createElement('div');
    geneStack.id = 'geneStack';
    geneStack.style.cssText = `
    position:absolute;
    right:14px;
    top:${ui.offsetTop + ui.offsetHeight + 10}px; 
    width:260px;
    max-height:40vh;
    overflow:auto;
    padding:10px 12px;
    background:rgba(20,20,28,.65);
    backdrop-filter:blur(4px);
    border:1px solid #333;
    border-radius:12px;
    color:#ddd;
    font:25px system-ui;
    line-height:1.4;
    z-index:999;
    display:none;
    `;
    document.body.appendChild(geneStack);

  const thrEl = ui.querySelector('#thr');
  const thrValEl = ui.querySelector('#thrVal');
  const showOnlyChk = ui.querySelector('#showOnlyChk');

  thrEl.oninput = () => {
    const thr = currentMin + parseFloat(thrEl.value) * (currentMax - currentMin);
    thrValEl.textContent = isFinite(thr) ? thr.toFixed(3) : '—';

    if (showOnlyChk.checked) {
      showOnlyExpressing(thr);
    } else {
      highlightPoints.visible = false;
      mainPoints.visible = true;
      //recolorByExpression(thr);
    }
    updateHighlight(thr); // soft glow for top cells either way
  };

  showOnlyChk.onchange = () => {
    const thr = currentMin + parseFloat(thrEl.value) * (currentMax - currentMin);
    if (showOnlyChk.checked) {
      showOnlyExpressing(thr);
    } else {
      showAllCells();
      //recolorByExpression(thr);
    }
  };
  
}



function showGeneStack(setName) {
  const genes = GENE_SETS[setName] || [];
  const stack = document.getElementById('geneStack');
  if (!stack) return;

  stack.innerHTML = genes.map(g => {
    const isMagic = MAGIC_GENES.has(g);
    const style = isMagic 
      ? 'color:#00e5ff;font-weight:700;text-shadow:0 0 6px rgba(0,229,255,.6);'
      : '';
    const prefix = isMagic ? '✨ ' : '';
    return `<div style="margin:2px 0;${style}">${prefix}${g}</div>`;
  }).join('');

  stack.style.display = 'block';
}

function hideGeneStack() {
  const stack = document.getElementById('geneStack');
  if (stack) stack.style.display = 'none';
}

function updateThrUI(thr){
  const uiVal = document.getElementById('thrVal');
  const uiSlider = document.getElementById('thr');
  if (!uiVal || !uiSlider) return;
  const span = (currentMax - currentMin) || 1;
  const t = (thr - currentMin) / span;
  uiSlider.value = String(Math.max(0, Math.min(1, t)));
  uiVal.textContent = isFinite(thr) ? thr.toFixed(3) : '—';
}

// === ANIMATION LOOP ===
function animate() {
  const ROT_SPEED_MAIN = 0.001;
  const ROT_SPEED_OVERLAY = 0.001;

  requestAnimationFrame(animate);

  // thumbnails (no bloom)
  for (let i = 0; i < NUM_THUMBS; i++) {
    const r = getThumbRect(i);
    const x = r.x + PANE_INNER_PAD;
    const y = r.y + PANE_INNER_PAD;
    const w = Math.max(1, r.w - PANE_INNER_PAD * 2);
    const h = Math.max(1, r.h - PANE_INNER_PAD * 2);
    renderer.setViewport(x, y, w, h);
    renderer.setScissor(x, y, w, h);

    if (thumbs[i].points) thumbs[i].points.rotation.y += 0.0015;
    renderer.render(thumbs[i].scene, thumbs[i].camera);
  }

  // main (with bloom)
  const mr = getMainRect();
  renderer.setViewport(mr.x, mr.y, mr.w, mr.h);
  renderer.setScissor(mr.x, mr.y, mr.w, mr.h);

  if (mainPoints) mainPoints.rotation.y += ROT_SPEED_MAIN;
  if (highlightPoints && highlightPoints.visible) { highlightPoints.rotation.y += ROT_SPEED_OVERLAY; }
  controls.update();
  composer.render();

  layoutLabels();
}
animate();
