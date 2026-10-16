// Shared by spacerocks' live pages (served as retro.js): formatting, and pixel drawing on small
// canvases with an NES-like palette and a 3x5 bitmap font.

const $ = id => document.getElementById(id);
const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

function calendar(jd) {
  const d = new Date((jd - 2440587.5) * 86400000);
  if (isNaN(d)) return "";
  const y = d.getUTCFullYear();
  if (y < 0 || y > 9999) return `year ${y.toLocaleString()}`;
  return d.toISOString().replace("T", " ").slice(0, 19);
}
function epochText(jd) { return `JD ${jd.toFixed(5)} · ${calendar(jd)}`; }
function duration(s) {
  if (!isFinite(s) || s < 0) return "–";
  if (s < 60) return `${s.toFixed(1)} s`;
  const h = Math.floor(s / 3600), m = Math.floor((s % 3600) / 60), sec = Math.floor(s % 60);
  return h ? `${h}h ${String(m).padStart(2, "0")}m` : `${m}m ${String(sec).padStart(2, "0")}s`;
}
function sci(x) { return x === 0 ? "0" : x.toExponential(2); }
function fmtDays(x) {
  const a = Math.abs(x);
  return a !== 0 && (a < 1e-3 || a >= 1e5) ? x.toExponential(3) : x.toPrecision(4);
}

// --- Pixel drawing -----------------------------------------------------------------------------
// A small NES-like palette and a 3x5 bitmap font, drawn one pixel at a time.
const PAL = {
  bg: "#0c0c24", grid: "#26265c", dim: "#5c5ca8", text: "#fcfcfc", star: "#3c3c7c",
  sun: "#f8b800", tp: "#7cd8fc", tracked: "#f83800", orbit: "#f87858",
  energy: "#58f898", series: ["#fca044", "#58f898", "#d878f8"],
  planets: ["#f87858", "#58d854", "#6888fc", "#f878f8", "#fca044", "#00e8d8", "#bcbcbc"],
};
const GLYPHS = {};
`0 111101101101111|1 010110010010111|2 111001111100111|3 111001111001111|4 101101111001001
5 111100111001111|6 111100111101111|7 111001001010010|8 111101111101111|9 111101111001111
A 010101111101101|B 110101110101110|C 011100100100011|D 110101101101110|E 111100110100111
F 111100110100100|G 011100101101011|H 101101111101101|I 111010010010111|J 001001001101010
K 101101110101101|L 100100100100111|M 101111111101101|N 110101101101101|O 010101101101010
P 110101110100100|Q 010101101110011|R 110101110101101|S 011100010001110|T 111010010010010
U 101101101101111|V 101101101101010|W 101101111111101|X 101101010101101|Y 101101010010010
Z 111001010100111|. 000000000000010|- 000000111000000|+ 000010111010000|: 000010000010000
/ 001001010100100|± 010111010000111|( 010100100100010|) 010001001001010|% 101001010100101
_ 000000000000111|' 010010000000000|= 000111000111000|° 010101010000000`.split(/[|\n]/).forEach(g => { GLYPHS[g[0]] = g.slice(2); });

function screen(id) {
  const c = $(id), ctx = c.getContext("2d");
  ctx.imageSmoothingEnabled = false;
  ctx.fillStyle = PAL.bg; ctx.fillRect(0, 0, c.width, c.height);
  return [ctx, c.width, c.height];
}
function px(ctx, x, y, color, size = 1) {
  ctx.fillStyle = color; ctx.fillRect(Math.round(x), Math.round(y), size, size);
}
function text(ctx, str, x, y, color = PAL.text) {
  ctx.fillStyle = color;
  x = Math.round(x); y = Math.round(y);
  for (const ch of String(str).toUpperCase()) {
    const g = GLYPHS[ch];
    if (g) for (let i = 0; i < 15; i++) if (g[i] === "1") ctx.fillRect(x + i % 3, y + Math.floor(i / 3), 1, 1);
    x += 4;
  }
}
const textWidth = str => String(str).length * 4 - 1;
function line(ctx, x0, y0, x1, y1, color) {
  // Bresenham, so lines stay one crisp pixel wide.
  x0 = Math.round(x0); y0 = Math.round(y0); x1 = Math.round(x1); y1 = Math.round(y1);
  const dx = Math.abs(x1 - x0), dy = -Math.abs(y1 - y0), sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
  let err = dx + dy;
  ctx.fillStyle = color;
  for (let n = 0; n < 2000; n++) {
    ctx.fillRect(x0, y0, 1, 1);
    if (x0 === x1 && y0 === y1) break;
    const e2 = 2 * err;
    if (e2 >= dy) { err += dy; x0 += sx; }
    if (e2 <= dx) { err += dx; y0 += sy; }
  }
}
function dotted(ctx, x0, y0, x1, y1, color) {
  ctx.fillStyle = color;
  if (y0 === y1) for (let x = x0; x <= x1; x += 2) ctx.fillRect(x, y0, 1, 1);
  else for (let y = y0; y <= y1; y += 2) ctx.fillRect(x0, y, 1, 1);
}
/// Short numbers for 4-pixel-wide characters: 3 significant figures, E notation when far from 1.
function short(x) {
  if (!isFinite(x)) return "-";
  const a = Math.abs(x);
  if (a !== 0 && (a < 1e-2 || a >= 1e5)) return x.toExponential(1).replace("e+", "E").replace("e", "E");
  return String(+x.toPrecision(3));
}
// A fixed starfield (a seeded generator, so it doesn't twinkle between frames).
const STARS = (() => { let r = 7; const n = () => (r = (r * 16807) % 2147483647) / 2147483647; return Array.from({ length: 70 }, () => [n(), n(), n() < 0.2]); })();
