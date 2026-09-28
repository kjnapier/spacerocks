//! Line-oriented driver used to cross-check the SPICE implementation against CSPICE.
//!
//! Commands (one per line on stdin):
//!   load <path>
//!   spkgeo <target> <observer> <et>
//!   state <target> <observer> <frame> <et>
//!   sxform <from> <to> <et>
//!   bary <jd_tdb> <id> <id> ...
//! Each command prints one line: `ok v1 v2 ...` or `err <message>`.

use spacerocks::spice::SpiceKernel;
use std::io::{self, BufRead, Write};

fn main() {
    let mut k = SpiceKernel::new();
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let p: Vec<&str> = line.split_whitespace().collect();
        if p.is_empty() {
            continue;
        }
        let res: Result<Vec<f64>, String> = match p[0] {
            "load" => k.load(p[1]).map(|_| vec![]).map_err(|e| e.to_string()),
            "unload" => Ok(vec![k.unload(p[1]) as i32 as f64]),
            "spkgeo" => k
                .spkgeo(p[1].parse().unwrap(), p[2].parse().unwrap(), p[3].parse().unwrap())
                .map(|s| s.to_vec())
                .map_err(|e| e.to_string()),
            "state" => k
                .state(p[1].parse().unwrap(), p[2].parse().unwrap(), p[3], p[4].parse().unwrap())
                .map(|s| s.to_vec())
                .map_err(|e| e.to_string()),
            "sxform" => k
                .sxform(p[1], p[2], p[3].parse().unwrap())
                .map(|x| {
                    let mut v = Vec::new();
                    for r in x.rotation.iter() {
                        v.extend_from_slice(r);
                    }
                    for r in x.rate.iter() {
                        v.extend_from_slice(r);
                    }
                    v
                })
                .map_err(|e| e.to_string()),
            "bary" => {
                let jd: f64 = p[1].parse().unwrap();
                let ids: Vec<i32> = p[2..].iter().map(|s| s.parse().unwrap()).collect();
                let mut o = vec![[0.0; 6]; ids.len()];
                k.barycentric_states_au(&ids, jd, &mut o)
                    .map(|_| o.iter().flatten().copied().collect())
                    .map_err(|e| e.to_string())
            }
            "bodyid" => k.body_id(&p[1..].join(" ")).map(|v| vec![v as f64]).ok_or("unknown".to_string()),
            other => Err(format!("unknown command {}", other)),
        };
        match res {
            Ok(v) => {
                let s: Vec<String> = v.iter().map(|x| format!("{:.17e}", x)).collect();
                writeln!(out, "ok {}", s.join(" ")).unwrap();
            }
            Err(e) => writeln!(out, "err {}", e.replace('\n', " ")).unwrap(),
        }
    }
}
