//! Text and notebook (HTML) views of what a `SpiceKernel` has loaded.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::sync::atomic::{AtomicUsize, Ordering};

use spacerocks::spice::{jd_from_et, KernelKind, KernelSummary, SegmentGroup, SpiceKernel};

// ------------------------------------------------------------------------------------------
// Small formatting helpers
// ------------------------------------------------------------------------------------------

/// Calendar date (TDB) of an ET, as (year, month, day). Julian calendar before 1582-10-15,
/// Gregorian after, as SPICE does by default.
fn ymd(et: f64) -> (i64, u32, u32) {
    let jd = jd_from_et(et) + 0.5;
    let z = jd.floor() as i64;
    let a = if z < 2_299_161 {
        z
    } else {
        let alpha = ((z as f64 - 1_867_216.25) / 36_524.25).floor() as i64;
        z + 1 + alpha - alpha.div_euclid(4)
    };
    let b = a + 1524;
    let c = ((b as f64 - 122.1) / 365.25).floor() as i64;
    let d = (365.25 * c as f64).floor() as i64;
    let e = ((b - d) as f64 / 30.6001).floor() as i64;
    let day = (b - d - (30.6001 * e as f64).floor() as i64) as u32;
    let month = if e < 14 { e - 1 } else { e - 13 } as u32;
    let year = if month > 2 { c - 4716 } else { c - 4715 };
    (year, month, day)
}

pub(crate) fn date(et: f64) -> String {
    let (y, m, d) = ymd(et);
    format!("{y:04}-{m:02}-{d:02}")
}

/// Approximate ET of 1 January of `year` (good to a day; used only for axis ticks).
fn et_of_year(year: i64) -> f64 {
    ((year - 2000) as f64 * 365.2425 - 0.5) * 86400.0
}

fn size(bytes: Option<u64>) -> String {
    match bytes {
        None => "—".into(),
        Some(b) if b < 1024 => format!("{b} B"),
        Some(b) if b < 1024 * 1024 => format!("{:.0} KB", b as f64 / 1024.0),
        Some(b) if b < 1024 * 1024 * 1024 => format!("{:.1} MB", b as f64 / 1048576.0),
        Some(b) => format!("{:.2} GB", b as f64 / 1073741824.0),
    }
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

fn basename(path: &str) -> &str {
    path.rsplit(['/', '\\']).next().unwrap_or(path)
}

fn dirname(path: &str) -> &str {
    match path.rfind(['/', '\\']) {
        Some(i) => &path[..i + 1],
        None => "",
    }
}

/// Current ET, approximately (37 leap seconds; used only to mark "today").
fn now_et() -> Option<f64> {
    let unix = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).ok()?.as_secs_f64();
    Some(unix - 946_727_935.816 + 69.184)
}

fn body_label(k: &SpiceKernel, id: i32) -> String {
    match k.body_name(id) {
        Some(n) => format!("{n} ({id})"),
        None => id.to_string(),
    }
}

fn frame_label(k: &SpiceKernel, id: i32) -> String {
    k.frame_info(id).map(|f| f.name).unwrap_or_else(|| format!("frame {id}"))
}

/// Name of the frame whose orientation a binary PCK class ID describes (the class ID is not
/// the frame ID: ITRF93 is frame 13000 but PCK class 3000).
pub(crate) fn pck_frame_label(k: &SpiceKernel, class_id: i32) -> String {
    for name in ["ITRF93", "MOON_PA_DE440", "MOON_PA_DE421", "MOON_PA_DE418", "MOON_PA", "EARTH_FIXED"] {
        if let Ok(f) = k.frame(name) {
            if f.class == 2 && f.class_id == class_id {
                return f.name;
            }
        }
    }
    match k.frame_info(class_id) {
        Some(f) if f.class == 2 && f.class_id == class_id => f.name,
        _ => format!("PCK class {class_id}"),
    }
}

fn span(groups: &[SegmentGroup]) -> Option<(f64, f64)> {
    let start = groups.iter().flat_map(|g| g.coverage.iter()).map(|i| i.start).fold(f64::INFINITY, f64::min);
    let end = groups.iter().flat_map(|g| g.coverage.iter()).map(|i| i.end).fold(f64::NEG_INFINITY, f64::max);
    (start <= end).then_some((start, end))
}

fn types(groups: &[SegmentGroup]) -> String {
    let mut t: Vec<i32> = groups.iter().map(|g| g.data_type).collect();
    t.sort_unstable();
    t.dedup();
    let t: Vec<String> = t.iter().map(|x| x.to_string()).collect();
    format!("type{} {}", if t.len() > 1 { "s" } else { "" }, t.join(", "))
}

fn distinct_bodies(groups: &[SegmentGroup]) -> Vec<i32> {
    let mut v: Vec<i32> = groups.iter().map(|g| g.body).collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// One-line description of a file's contents.
pub(crate) fn describe(k: &SpiceKernel, s: &KernelSummary) -> String {
    match s.kind {
        KernelKind::Spk => {
            let bodies = distinct_bodies(&s.groups);
            if s.groups.len() == 1 {
                let g = &s.groups[0];
                format!(
                    "{} relative to {}, {}, {} segment{}",
                    body_label(k, g.body),
                    body_label(k, g.center),
                    types(&s.groups),
                    g.n_segments,
                    if g.n_segments == 1 { "" } else { "s" }
                )
            } else {
                format!("{} bodies, {}", bodies.len(), types(&s.groups))
            }
        }
        KernelKind::Pck => {
            let frames: Vec<String> = distinct_bodies(&s.groups).iter().map(|&c| pck_frame_label(k, c)).collect();
            let refs: Vec<String> = {
                let mut r: Vec<i32> = s.groups.iter().map(|g| g.frame).collect();
                r.sort_unstable();
                r.dedup();
                r.iter().map(|&f| frame_label(k, f)).collect()
            };
            let n: usize = s.groups.iter().map(|g| g.n_segments).sum();
            format!("{} orientation relative to {}, {}, {} segments", frames.join(", "), refs.join(", "), types(&s.groups), n)
        }
        KernelKind::Text | KernelKind::Meta => {
            let v = &s.variables;
            let mut parts = Vec::new();
            if v.iter().any(|n| n.starts_with("DELTET/")) {
                parts.push("leapseconds".to_string());
            }
            let mut body_ids: Vec<&str> = v
                .iter()
                .filter_map(|n| n.strip_prefix("BODY"))
                .filter_map(|r| r.split('_').next())
                .filter(|id| !id.is_empty() && id.trim_start_matches('-').chars().all(|c| c.is_ascii_digit()))
                .collect();
            body_ids.sort_unstable();
            body_ids.dedup();
            if !body_ids.is_empty() {
                parts.push(format!("constants for {} bod{}", body_ids.len(), if body_ids.len() == 1 { "y" } else { "ies" }));
            }
            let n_frames = v.iter().filter(|n| n.starts_with("FRAME_") && n.ends_with("_NAME")).count();
            if n_frames > 0 {
                parts.push(format!("{n_frames} frame{}", if n_frames == 1 { "" } else { "s" }));
            }
            if v.iter().any(|n| n == "NAIF_BODY_NAME") {
                parts.push("body names".into());
            }
            if s.kind == KernelKind::Meta {
                parts.insert(0, format!("loads {} file{}", s.children.len(), if s.children.len() == 1 { "" } else { "s" }));
            }
            let vars = format!("{} variable{}", v.len(), if v.len() == 1 { "" } else { "s" });
            if parts.is_empty() {
                vars
            } else {
                format!("{} ({vars})", parts.join(", "))
            }
        }
    }
}

fn coverage_text(s: &KernelSummary) -> String {
    match span(&s.groups) {
        Some((a, b)) => format!("{} → {}", date(a), date(b)),
        None => String::new(),
    }
}

// ------------------------------------------------------------------------------------------
// Plain text
// ------------------------------------------------------------------------------------------

pub(crate) fn text(k: &SpiceKernel) -> String {
    let files = k.summary();
    if files.is_empty() {
        return "SpiceKernel(no kernels loaded)".to_string();
    }
    let mut s = format!("SpiceKernel ({} file{}, highest priority last):\n", files.len(), if files.len() == 1 { "" } else { "s" });
    for f in &files {
        let cov = coverage_text(f);
        let _ = write!(s, "  [{}] {}\n        {}", f.kind, f.path, describe(k, f));
        if !cov.is_empty() {
            let _ = write!(s, "; {cov}");
        }
        s.push('\n');
    }
    s
}

// ------------------------------------------------------------------------------------------
// HTML (Jupyter `_repr_html_`)
// ------------------------------------------------------------------------------------------

/// Categorical file colours (light, dark); fixed order, never cycled — files past the
/// eighth share a neutral grey.
const PALETTE: [(&str, &str); 8] = [
    ("#2a78d6", "#3987e5"),
    ("#eb6834", "#d95926"),
    ("#1baf7a", "#199e70"),
    ("#eda100", "#c98500"),
    ("#e87ba4", "#d55181"),
    ("#008300", "#008300"),
    ("#4a3aa7", "#9085e9"),
    ("#e34948", "#e66767"),
];
const OTHER: (&str, &str) = ("#8a8984", "#8a8984");

static COUNTER: AtomicUsize = AtomicUsize::new(0);

fn style(root: &str) -> String {
    let mut vars_light = String::new();
    let mut vars_dark = String::new();
    for (i, (l, d)) in PALETTE.iter().chain(std::iter::once(&OTHER)).enumerate() {
        let _ = write!(vars_light, "--c{i}:{l};");
        let _ = write!(vars_dark, "--c{i}:{d};");
    }
    let dark = format!(
        "{vars_dark}--ink:#f2f1ed;--ink2:#c3c2b7;--ink3:#8f8e87;--rule:rgba(255,255,255,.14);--band:rgba(255,255,255,.04);--gap:#1e1e1d;"
    );
    format!(
        r#"<style>
.{root}{{{vars_light}--ink:#0b0b0b;--ink2:#52514e;--ink3:#8a8984;--rule:rgba(0,0,0,.12);--band:rgba(0,0,0,.03);--gap:#ffffff;
  font-family:var(--jp-ui-font-family,-apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif);font-size:12.5px;line-height:1.45;color:var(--ink);max-width:980px}}
@media (prefers-color-scheme: dark){{.{root}{{{dark}}}}}
body[data-jp-theme-light="false"] .{root},body.vscode-dark .{root},body.vscode-high-contrast .{root},:root[data-theme="dark"] .{root}{{{dark}}}
body[data-jp-theme-light="true"] .{root},body.vscode-light .{root},:root[data-theme="light"] .{root}{{{vars_light}--ink:#0b0b0b;--ink2:#52514e;--ink3:#8a8984;--rule:rgba(0,0,0,.12);--band:rgba(0,0,0,.03);--gap:#ffffff}}
.{root} .hd{{display:flex;gap:14px;align-items:baseline;flex-wrap:wrap;margin:2px 0 8px}}
.{root} .hd b{{font-size:14px;font-weight:600}}
.{root} .hd span{{color:var(--ink2)}}
.{root} table{{border-collapse:collapse;width:100%;margin:0 0 10px;background:transparent}}
.{root} th{{text-align:left;font-weight:500;color:var(--ink3);font-size:11px;text-transform:uppercase;letter-spacing:.04em;padding:3px 10px 5px 0;border-bottom:1px solid var(--rule);background:transparent}}
.{root} td{{text-align:left;vertical-align:top;padding:5px 10px 5px 0;border-bottom:1px solid var(--rule);background:transparent;color:var(--ink)}}
.{root} tr:hover td{{background:var(--band)}}
.{root} td.n,.{root} th.n{{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}}
.{root} td.dim{{color:var(--ink2)}}
.{root} .sw{{display:inline-block;width:10px;height:10px;border-radius:3px;vertical-align:-1px}}
.{root} .sw.none{{border:1px solid var(--rule);width:8px;height:8px}}
.{root} .kind{{font-size:10.5px;font-weight:600;letter-spacing:.03em;color:var(--ink2);border:1px solid var(--rule);border-radius:4px;padding:0 5px;white-space:nowrap}}
.{root} .file{{font-family:var(--jp-code-font-family,ui-monospace,SFMono-Regular,Menlo,monospace);font-size:12px}}
.{root} .path{{color:var(--ink3);font-size:11px;word-break:break-all}}
.{root} .ttl{{font-size:11px;text-transform:uppercase;letter-spacing:.04em;color:var(--ink3);margin:6px 0 2px}}
.{root} svg text{{font-family:inherit;fill:var(--ink2);font-size:11px}}
.{root} svg .lbl{{fill:var(--ink)}}
.{root} svg .tick{{fill:var(--ink3);font-size:10px}}
.{root} svg .ax{{stroke:var(--rule)}}
.{root} svg g.bar:hover{{opacity:.8}}
.{root} .tl>input{{position:absolute;opacity:0;pointer-events:none}}
.{root} .tl>label{{display:inline-block;cursor:pointer;font-size:11.5px;color:var(--ink2);border:1px solid var(--rule);padding:1px 9px;margin:2px 0 6px}}
.{root} .tl>label:first-of-type{{border-radius:6px 0 0 6px}}
.{root} .tl>label:last-of-type{{border-radius:0 6px 6px 0;border-left:none}}
.{root} .tl>input:checked+label{{color:var(--ink);background:var(--band);font-weight:600}}
.{root} .tl>input:focus-visible+label{{outline:2px solid var(--c0)}}
.{root} .tl .hint{{color:var(--ink3);font-size:11px;margin-left:10px}}
.{root} .tl .v-a{{display:none}}
.{root} .tl>input:last-of-type:checked~.v-f{{display:none}}
.{root} .tl>input:last-of-type:checked~.v-a{{display:block}}
.{root} .tl>input:last-of-type:checked~.hint{{display:none}}
.{root} details{{margin-top:6px}}
.{root} summary{{cursor:pointer;color:var(--ink2)}}
.{root} details table td,.{root} details table th{{padding-top:3px;padding-bottom:3px}}
</style>"#
    )
}

struct Row {
    label: String,
    /// (file index, group) — lowest priority first so higher-priority bars draw on top.
    bars: Vec<(usize, SegmentGroup, String)>,
}

fn nice_step(years: f64) -> i64 {
    for s in [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 5000] {
        if years / s as f64 <= 8.0 {
            return s;
        }
    }
    10000
}


/// Days from 1970-01-01 to a proleptic Gregorian date.
fn days_from_civil(y: i64, m: u32, d: u32) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let m = m as i64;
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + d as i64 - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// ET at 00:00 of a calendar date (TDB, to within leap-second-sized error; used for ticks).
fn et_of_date(y: i64, m: u32, d: u32) -> f64 {
    (days_from_civil(y, m, d) - days_from_civil(2000, 1, 1)) as f64 * 86400.0 - 43200.0
}

/// Axis ticks (ET, label) for the range [t0, t1]: years, months or days, at most ~9.
fn ticks(t0: f64, t1: f64) -> Vec<(f64, String)> {
    let days = (t1 - t0) / 86400.0;
    let mut out = Vec::new();
    if days > 3.0 * 365.25 {
        let (y0, y1) = (ymd(t0).0, ymd(t1).0);
        let step = nice_step((y1 - y0).max(1) as f64);
        let mut y = y0.div_euclid(step) * step;
        while y <= y1 + step {
            let e = et_of_year(y);
            if e >= t0 && e <= t1 {
                out.push((e, y.to_string()));
            }
            y += step;
        }
    } else if days > 60.0 {
        let step = [1u32, 2, 3, 6, 12].into_iter().find(|s| days / 30.44 / *s as f64 <= 8.0).unwrap_or(12);
        let (mut y, mut m, _) = ymd(t0);
        m = ((m - 1) / step) * step + 1;
        loop {
            let e = et_of_date(y, m, 1);
            if e > t1 {
                break;
            }
            if e >= t0 {
                out.push((e, format!("{y}-{m:02}")));
            }
            m += step;
            if m > 12 {
                m -= 12;
                y += 1;
            }
        }
    } else {
        let step = [1i64, 2, 5, 7, 10, 14].into_iter().find(|s| days / *s as f64 <= 8.0).unwrap_or(14);
        let (y, m, d) = ymd(t0);
        let mut day = days_from_civil(y, m, d);
        while (day as f64 - days_from_civil(2000, 1, 1) as f64) * 86400.0 - 43200.0 <= t1 {
            let e = (day - days_from_civil(2000, 1, 1)) as f64 * 86400.0 - 43200.0;
            if e >= t0 {
                out.push((e, date(e)));
            }
            day += step;
        }
    }
    out
}

/// When a few rows span most of the axis and squash the rest, the range to zoom to: the
/// extent of the short rows (and today, if it falls inside), padded. None if no zoom helps.
fn focus_range(spans: &[(f64, f64)], t0: f64, t1: f64) -> Option<(f64, f64)> {
    let full = t1 - t0;
    let short: Vec<&(f64, f64)> = spans.iter().filter(|s| s.1 - s.0 < 0.25 * full).collect();
    if short.is_empty() {
        return None;
    }
    let mut f0 = short.iter().map(|s| s.0).fold(f64::INFINITY, f64::min);
    let mut f1 = short.iter().map(|s| s.1).fold(f64::NEG_INFINITY, f64::max);
    if let Some(now) = now_et().filter(|&t| t > t0 && t < t1) {
        f0 = f0.min(now);
        f1 = f1.max(now);
    }
    let pad = ((f1 - f0) * 0.06).max(86400.0);
    let (f0, f1) = ((f0 - pad).max(t0), (f1 + pad).min(t1));
    ((f1 - f0) < 0.5 * full).then_some((f0, f1))
}

/// SVG timeline of `rows` over the axis [t0, t1]; bars outside are clipped and marked with
/// an arrow at the edge they continue past. Row labels on the right give the true extent.
fn timeline(rows: &[Row], slot: &[Option<usize>], t0: f64, t1: f64) -> String {
    let (w, lw, rw, rh, top) = (960.0, 230.0, 150.0, 20.0, 22.0);
    let pw = w - lw - rw - 16.0;
    let x = |et: f64| lw + (et.clamp(t0, t1) - t0) / (t1 - t0) * pw;
    let height = top + rows.len() as f64 * rh + 6.0;
    let mut h = String::new();
    let _ = write!(
        h,
        r#"<div style="overflow-x:auto"><svg viewBox="0 0 {w} {height}" width="100%" style="max-width:{w}px;min-width:720px;display:block" role="img" aria-label="Coverage timeline of loaded ephemerides and orientation data">"#
    );
    for (e, label) in ticks(t0, t1) {
        let xx = x(e);
        let _ = write!(
            h,
            r#"<line class="ax" x1="{xx:.1}" x2="{xx:.1}" y1="{}" y2="{:.1}" stroke-width="1"/><text class="tick" x="{xx:.1}" y="12" text-anchor="middle">{label}</text>"#,
            top - 4.0,
            height - 4.0
        );
    }
    for (r, row) in rows.iter().enumerate() {
        let yc = top + r as f64 * rh + rh / 2.0;
        let _ = write!(h, r#"<text class="lbl" x="{:.1}" y="{:.1}" text-anchor="end" dominant-baseline="middle">{}</text>"#, lw - 12.0, yc, esc(&row.label));
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for (fi, g, tip) in &row.bars {
            let c = slot[*fi].unwrap_or(PALETTE.len());
            for iv in &g.coverage {
                lo = lo.min(iv.start);
                hi = hi.max(iv.end);
                if iv.end < t0 || iv.start > t1 {
                    continue;
                }
                let (cl, cr) = (iv.start < t0, iv.end > t1);
                let xa = x(iv.start);
                let bw = (x(iv.end) - xa).max(3.0);
                let y = yc - 5.0;
                let title = format!("<title>{}&#10;{} → {}</title>", esc(tip), date(iv.start), date(iv.end));
                let _ = write!(
                    h,
                    r#"<g class="bar"><rect x="{xa:.1}" y="{y:.1}" width="{bw:.1}" height="10" rx="3" fill="var(--c{c})" stroke="var(--gap)" stroke-width="1.5" paint-order="stroke"/>"#
                );
                // Square off clipped ends and point an arrow past the axis.
                if cl {
                    let _ = write!(
                        h,
                        r#"<rect x="{xa:.1}" y="{y:.1}" width="4" height="10" fill="var(--c{c})"/><path d="M{:.1} {:.1}l-5 5l5 5z" fill="var(--c{c})"/>"#,
                        xa - 2.0,
                        y
                    );
                }
                if cr {
                    let xe = xa + bw;
                    let _ = write!(
                        h,
                        r#"<rect x="{:.1}" y="{y:.1}" width="4" height="10" fill="var(--c{c})"/><path d="M{:.1} {:.1}l5 5l-5 5z" fill="var(--c{c})"/>"#,
                        xe - 4.0,
                        xe + 2.0,
                        y
                    );
                }
                let _ = write!(h, "{title}</g>");
            }
        }
        let _ = write!(h, r#"<text x="{:.1}" y="{:.1}" dominant-baseline="middle">{} → {}</text>"#, w - rw, yc, date(lo), date(hi));
    }
    if let Some(now) = now_et().filter(|&t| t > t0 && t < t1) {
        let xx = x(now);
        let _ = write!(
            h,
            r#"<line x1="{xx:.1}" x2="{xx:.1}" y1="{}" y2="{:.1}" stroke="var(--ink2)" stroke-width="1" stroke-dasharray="3 3"/><text x="{:.1}" y="{:.1}" style="font-size:10px">today</text>"#,
            top - 4.0,
            height - 4.0,
            xx + 3.0,
            top + 2.0
        );
    }
    h.push_str("</svg></div>");
    h
}

const MAX_ROWS: usize = 60;

pub(crate) fn html(k: &SpiceKernel) -> String {
    let root = format!("srk{}", COUNTER.fetch_add(1, Ordering::Relaxed));
    let files = k.summary();
    let mut h = String::new();
    h.push_str(&style(&root));
    let _ = write!(h, r#"<div class="{root}">"#);

    if files.is_empty() {
        h.push_str(r#"<div class="hd"><b>SpiceKernel</b><span>no kernels loaded</span></div></div>"#);
        return h;
    }

    // Colour slot per binary file, in load order; text kernels get none.
    let mut slot = vec![None; files.len()];
    let mut next = 0usize;
    for (i, f) in files.iter().enumerate() {
        if matches!(f.kind, KernelKind::Spk | KernelKind::Pck) {
            slot[i] = Some(next.min(PALETTE.len()));
            next += 1;
        }
    }

    // Header.
    let n_bodies = k.spk_bodies().len();
    let n_pool = k.pool().len();
    let n_frames: usize = {
        let mut ids: Vec<i32> = files.iter().filter(|f| f.kind == KernelKind::Pck).flat_map(|f| f.groups.iter().map(|g| g.body)).collect();
        ids.sort_unstable();
        ids.dedup();
        ids.len()
    };
    let _ = write!(
        h,
        r#"<div class="hd"><b>SpiceKernel</b><span>{} file{}</span><span>{} bod{} with ephemerides</span>"#,
        files.len(),
        if files.len() == 1 { "" } else { "s" },
        n_bodies,
        if n_bodies == 1 { "y" } else { "ies" }
    );
    if n_frames > 0 {
        let _ = write!(h, "<span>{n_frames} body-fixed frame{}</span>", if n_frames == 1 { "" } else { "s" });
    }
    if n_pool > 0 {
        let _ = write!(h, "<span>{n_pool} pool variable{}</span>", if n_pool == 1 { "" } else { "s" });
    }
    h.push_str("</div>");

    // File table, highest priority first.
    h.push_str(r#"<table><thead><tr><th class="n">#</th><th></th><th>Type</th><th>File</th><th>Contents</th><th>Coverage (TDB)</th><th class="n">Size</th></tr></thead><tbody>"#);
    for (rank, (i, f)) in files.iter().enumerate().rev().enumerate() {
        let sw = match slot[i] {
            Some(c) => format!(r#"<span class="sw" style="background:var(--c{c})"></span>"#),
            None => r#"<span class="sw none"></span>"#.to_string(),
        };
        let kind = match f.kind {
            KernelKind::Spk => "SPK",
            KernelKind::Pck => "PCK",
            KernelKind::Text => "TEXT",
            KernelKind::Meta => "META",
        };
        let _ = write!(
            h,
            r#"<tr><td class="n dim">{}</td><td>{sw}</td><td><span class="kind">{kind}</span></td><td><div class="file">{}</div><div class="path" title="{}">{}</div></td><td>{}</td><td class="dim" style="white-space:nowrap">{}</td><td class="n dim">{}</td></tr>"#,
            rank + 1,
            esc(basename(&f.path)),
            esc(&f.path),
            esc(dirname(&f.path)),
            esc(&describe(k, f)),
            coverage_text(f),
            size(f.size_bytes)
        );
    }
    h.push_str("</tbody></table>");

    // Timeline rows: one per SPK body, then one per PCK frame.
    let mut body_rows: BTreeMap<(bool, i32), Row> = BTreeMap::new();
    let mut frame_rows: BTreeMap<i32, Row> = BTreeMap::new();
    for (i, f) in files.iter().enumerate() {
        for g in &f.groups {
            let tip = match f.kind {
                KernelKind::Spk => format!(
                    "{} relative to {} · {} · type {} · {} segment{}",
                    body_label(k, g.body),
                    body_label(k, g.center),
                    basename(&f.path),
                    g.data_type,
                    g.n_segments,
                    if g.n_segments == 1 { "" } else { "s" }
                ),
                _ => format!(
                    "{} relative to {} · {} · type {} · {} segment{}",
                    pck_frame_label(k, g.body),
                    frame_label(k, g.frame),
                    basename(&f.path),
                    g.data_type,
                    g.n_segments,
                    if g.n_segments == 1 { "" } else { "s" }
                ),
            };
            match f.kind {
                KernelKind::Spk => body_rows
                    .entry((g.body < 0, if g.body < 0 { -g.body } else { g.body }))
                    .or_insert_with(|| Row { label: body_label(k, g.body), bars: Vec::new() })
                    .bars
                    .push((i, g.clone(), tip)),
                KernelKind::Pck => frame_rows
                    .entry(g.body)
                    .or_insert_with(|| Row { label: format!("{} orientation", pck_frame_label(k, g.body)), bars: Vec::new() })
                    .bars
                    .push((i, g.clone(), tip)),
                _ => {}
            }
        }
    }
    let mut rows: Vec<Row> = body_rows.into_values().chain(frame_rows.into_values()).collect();

    if !rows.is_empty() {
        let hidden = rows.len().saturating_sub(MAX_ROWS);
        rows.truncate(MAX_ROWS);

        // Full extent of everything plotted.
        let row_span = |r: &Row| {
            let lo = r.bars.iter().flat_map(|b| b.1.coverage.iter()).map(|i| i.start).fold(f64::INFINITY, f64::min);
            let hi = r.bars.iter().flat_map(|b| b.1.coverage.iter()).map(|i| i.end).fold(f64::NEG_INFINITY, f64::max);
            (lo, hi)
        };
        let spans: Vec<(f64, f64)> = rows.iter().map(row_span).collect();
        let t0 = spans.iter().map(|s| s.0).fold(f64::INFINITY, f64::min);
        let t1 = spans.iter().map(|s| s.1).fold(f64::NEG_INFINITY, f64::max);
        let (t0, t1) = if t1 > t0 { (t0, t1) } else { (t0 - 86400.0, t1 + 86400.0) };

        h.push_str(r#"<div class="ttl">Coverage</div>"#);
        match focus_range(&spans, t0, t1) {
            Some((f0, f1)) => {
                let (fy0, fy1, ay0, ay1) = (ymd(f0).0, ymd(f1).0, ymd(t0).0, ymd(t1).0);
                let _ = write!(
                    h,
                    r#"<div class="tl"><input type="radio" name="{root}-tl" id="{root}-f" checked><label for="{root}-f">Zoomed {fy0}–{fy1}</label><input type="radio" name="{root}-tl" id="{root}-a"><label for="{root}-a">Full range {ay0}–{ay1}</label><span class="hint">Bars with arrows continue past the axis.</span><div class="v-f">{}</div><div class="v-a">{}</div></div>"#,
                    timeline(&rows, &slot, f0, f1),
                    timeline(&rows, &slot, t0, t1)
                );
            }
            None => h.push_str(&timeline(&rows, &slot, t0, t1)),
        }
        if hidden > 0 {
            let _ = write!(h, r#"<div class="path">…and {hidden} more (see the segment table below)</div>"#);
        }

        // Segment table.
        h.push_str(r#"<details><summary>Segment table</summary><table><thead><tr><th></th><th>Body / frame</th><th>Relative to</th><th>Frame</th><th class="n">Type</th><th class="n">Segments</th><th>Coverage (TDB)</th><th>File</th></tr></thead><tbody>"#);
        for (i, f) in files.iter().enumerate().rev() {
            for g in &f.groups {
                let c = slot[i].unwrap_or(PALETTE.len());
                let (what, rel, frame) = match f.kind {
                    KernelKind::Spk => (body_label(k, g.body), body_label(k, g.center), frame_label(k, g.frame)),
                    _ => (format!("{} orientation", pck_frame_label(k, g.body)), frame_label(k, g.frame), "—".to_string()),
                };
                let cov: Vec<String> = g.coverage.iter().map(|iv| format!("{} → {}", date(iv.start), date(iv.end))).collect();
                let _ = write!(
                    h,
                    r#"<tr><td><span class="sw" style="background:var(--c{c})"></span></td><td>{}</td><td class="dim">{}</td><td class="dim">{}</td><td class="n">{}</td><td class="n">{}</td><td class="dim" style="white-space:nowrap">{}</td><td class="file">{}</td></tr>"#,
                    esc(&what),
                    esc(&rel),
                    esc(&frame),
                    g.data_type,
                    g.n_segments,
                    cov.join("<br>"),
                    esc(basename(&f.path))
                );
            }
        }
        h.push_str("</tbody></table></details>");
    }

    h.push_str("</div>");
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tick_granularity() {
        let d = 86400.0;
        let y = ticks(et_of_date(1850, 1, 1), et_of_date(2150, 1, 1));
        assert!(y.len() >= 4 && y.len() <= 9 && y.iter().all(|t| t.1.len() == 4), "{y:?}");
        let m = ticks(et_of_date(2026, 1, 10), et_of_date(2026, 9, 20));
        assert_eq!(m.iter().map(|t| t.1.as_str()).collect::<Vec<_>>(), ["2026-03", "2026-05", "2026-07", "2026-09"]);
        let m = ticks(et_of_date(2026, 1, 10), et_of_date(2026, 6, 20));
        assert_eq!(m.iter().map(|t| t.1.as_str()).collect::<Vec<_>>(), ["2026-02", "2026-03", "2026-04", "2026-05", "2026-06"]);
        let days = ticks(et_of_date(2026, 9, 1) + 3600.0, et_of_date(2026, 9, 1) + 20.0 * d);
        assert_eq!(days.iter().map(|t| t.1.as_str()).collect::<Vec<_>>(), ["2026-09-06", "2026-09-11", "2026-09-16", "2026-09-21"]);
        assert!(days.len() <= 9);
        assert_eq!(date(et_of_date(2026, 9, 28)), "2026-09-28");
    }

    #[test]
    fn focus_only_when_it_helps() {
        let yr = 365.25 * 86400.0;
        let long = (-150.0 * yr, 150.0 * yr);
        // One short row among long ones: zoom around it.
        let (f0, f1) = focus_range(&[long, long, (22.0 * yr, 31.0 * yr)], long.0, long.1).unwrap();
        assert!(f0 < 22.0 * yr && f1 > 31.0 * yr && f1 - f0 < 0.5 * (long.1 - long.0));
        // All rows comparable: no zoom.
        assert!(focus_range(&[long, (-100.0 * yr, 140.0 * yr)], long.0, long.1).is_none());
    }
}
