//! Text kernel parsing (LSK, text PCK, FK, SCLK, meta-kernels) and the kernel variable pool.

use std::collections::HashMap;

use super::error::{Result, SpiceError};

/// A value in the kernel pool: either all numeric or all strings.
#[derive(Debug, Clone, PartialEq)]
pub enum PoolValue {
    Numeric(Vec<f64>),
    Strings(Vec<String>),
}

/// The kernel variable pool populated from text kernels.
#[derive(Debug, Clone, Default)]
pub struct KernelPool {
    vars: HashMap<String, PoolValue>,
}

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Name(String),
    Eq,
    PlusEq,
    LParen,
    RParen,
    Num(f64),
    Str(String),
}

impl KernelPool {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn get(&self, name: &str) -> Option<&PoolValue> {
        self.vars.get(name)
    }

    pub fn contains(&self, name: &str) -> bool {
        self.vars.contains_key(name)
    }

    pub fn names(&self) -> impl Iterator<Item = &String> {
        self.vars.keys()
    }

    pub fn len(&self) -> usize {
        self.vars.len()
    }

    pub fn is_empty(&self) -> bool {
        self.vars.is_empty()
    }

    /// Numeric values of a variable, if it exists and is numeric.
    pub fn get_f64s(&self, name: &str) -> Option<&[f64]> {
        match self.vars.get(name) {
            Some(PoolValue::Numeric(v)) => Some(v),
            _ => None,
        }
    }

    /// First numeric value of a variable.
    pub fn get_f64(&self, name: &str) -> Option<f64> {
        self.get_f64s(name).and_then(|v| v.first().copied())
    }

    /// First value of a numeric variable, rounded to an integer.
    pub fn get_i32(&self, name: &str) -> Option<i32> {
        self.get_f64(name).map(|v| v.round() as i32)
    }

    /// String values of a variable, if it exists and holds strings.
    pub fn get_strs(&self, name: &str) -> Option<&[String]> {
        match self.vars.get(name) {
            Some(PoolValue::Strings(v)) => Some(v),
            _ => None,
        }
    }

    pub fn get_str(&self, name: &str) -> Option<&str> {
        self.get_strs(name).and_then(|v| v.first()).map(|s| s.as_str())
    }

    /// String values with SPICE continuation handling: an element ending in `+` is joined
    /// with the following element (as done by `STPOOL` for meta-kernels).
    pub fn get_strs_joined(&self, name: &str) -> Option<Vec<String>> {
        let raw = self.get_strs(name)?;
        let mut out = Vec::new();
        let mut cur = String::new();
        for s in raw {
            let t = s.trim_end();
            if let Some(stripped) = t.strip_suffix('+') {
                cur.push_str(stripped);
            } else {
                cur.push_str(t);
                out.push(std::mem::take(&mut cur));
            }
        }
        if !cur.is_empty() {
            out.push(cur);
        }
        Some(out)
    }

    /// Set (replace) a numeric variable.
    pub fn set_f64s(&mut self, name: &str, values: Vec<f64>) {
        self.vars.insert(name.to_string(), PoolValue::Numeric(values));
    }

    /// Set (replace) a string variable.
    pub fn set_strs(&mut self, name: &str, values: Vec<String>) {
        self.vars.insert(name.to_string(), PoolValue::Strings(values));
    }

    /// Parse text kernel contents and add its assignments to the pool.
    /// Returns the names of the variables that were assigned.
    pub fn load_str(&mut self, text: &str, path: &str) -> Result<Vec<String>> {
        let toks = tokenize(text, path)?;
        let mut assigned = Vec::new();
        let mut i = 0;
        while i < toks.len() {
            let (line, name) = match &toks[i] {
                (l, Tok::Name(n)) => (*l, n.clone()),
                (l, t) => {
                    return Err(SpiceError::TextKernel {
                        path: path.to_string(),
                        line: *l,
                        reason: format!("expected a variable name, found {:?}", t),
                    })
                }
            };
            i += 1;
            let append = match toks.get(i) {
                Some((_, Tok::Eq)) => false,
                Some((_, Tok::PlusEq)) => true,
                _ => {
                    return Err(SpiceError::TextKernel {
                        path: path.to_string(),
                        line,
                        reason: format!("expected '=' or '+=' after '{}'", name),
                    })
                }
            };
            i += 1;
            let mut nums = Vec::new();
            let mut strs = Vec::new();
            let mut push = |t: &Tok, line: usize| -> Result<()> {
                match t {
                    Tok::Num(v) => nums.push(*v),
                    Tok::Str(s) => strs.push(s.clone()),
                    other => {
                        return Err(SpiceError::TextKernel {
                            path: path.to_string(),
                            line,
                            reason: format!("unexpected token {:?} in value of '{}'", other, name),
                        })
                    }
                }
                Ok(())
            };
            match toks.get(i) {
                Some((_, Tok::LParen)) => {
                    i += 1;
                    loop {
                        match toks.get(i) {
                            Some((_, Tok::RParen)) => {
                                i += 1;
                                break;
                            }
                            Some((l, t)) => {
                                push(t, *l)?;
                                i += 1;
                            }
                            None => {
                                return Err(SpiceError::TextKernel {
                                    path: path.to_string(),
                                    line,
                                    reason: format!("unterminated value list for '{}'", name),
                                })
                            }
                        }
                    }
                }
                Some((l, t)) => {
                    push(t, *l)?;
                    i += 1;
                }
                None => {
                    return Err(SpiceError::TextKernel {
                        path: path.to_string(),
                        line,
                        reason: format!("missing value for '{}'", name),
                    })
                }
            }
            if !nums.is_empty() && !strs.is_empty() {
                return Err(SpiceError::TextKernel {
                    path: path.to_string(),
                    line,
                    reason: format!("'{}' mixes numeric and string values", name),
                });
            }
            let value = if !strs.is_empty() {
                PoolValue::Strings(strs)
            } else {
                PoolValue::Numeric(nums)
            };
            match (append, self.vars.get_mut(&name), value) {
                (true, Some(PoolValue::Numeric(v)), PoolValue::Numeric(new)) => v.extend(new),
                (true, Some(PoolValue::Strings(v)), PoolValue::Strings(new)) => v.extend(new),
                (true, Some(_), _) => {
                    return Err(SpiceError::TextKernel {
                        path: path.to_string(),
                        line,
                        reason: format!("'+=' changes the type of '{}'", name),
                    })
                }
                (_, _, value) => {
                    self.vars.insert(name.clone(), value);
                }
            }
            assigned.push(name);
        }
        Ok(assigned)
    }
}

/// Split the data sections of a text kernel into tokens (with line numbers).
fn tokenize(text: &str, path: &str) -> Result<Vec<(usize, Tok)>> {
    let mut toks = Vec::new();
    let mut in_data = false;
    for (lineno, raw_line) in text.lines().enumerate() {
        let lineno = lineno + 1;
        let line = raw_line.replace('\t', " ").replace('\r', "");
        let trimmed = line.trim();
        // As in CSPICE, a marker counts only alone on its line: comment text such as
        // "\\begindata token. In order to ..." (pck00010.tpc) is not one.
        if trimmed == "\\begindata" {
            in_data = true;
            continue;
        }
        if trimmed == "\\begintext" {
            in_data = false;
            continue;
        }
        if !in_data {
            continue;
        }
        let chars: Vec<char> = line.chars().collect();
        let mut k = 0;
        while k < chars.len() {
            let c = chars[k];
            if c.is_whitespace() || c == ',' {
                k += 1;
                continue;
            }
            match c {
                '(' => {
                    toks.push((lineno, Tok::LParen));
                    k += 1;
                }
                ')' => {
                    toks.push((lineno, Tok::RParen));
                    k += 1;
                }
                '=' => {
                    toks.push((lineno, Tok::Eq));
                    k += 1;
                }
                '+' if k + 1 < chars.len() && chars[k + 1] == '=' => {
                    toks.push((lineno, Tok::PlusEq));
                    k += 2;
                }
                '\'' => {
                    let mut s = String::new();
                    k += 1;
                    let mut closed = false;
                    while k < chars.len() {
                        if chars[k] == '\'' {
                            if k + 1 < chars.len() && chars[k + 1] == '\'' {
                                s.push('\'');
                                k += 2;
                                continue;
                            }
                            closed = true;
                            k += 1;
                            break;
                        }
                        s.push(chars[k]);
                        k += 1;
                    }
                    if !closed {
                        return Err(SpiceError::TextKernel {
                            path: path.to_string(),
                            line: lineno,
                            reason: "unterminated string".into(),
                        });
                    }
                    toks.push((lineno, Tok::Str(s)));
                }
                '@' => {
                    let start = k + 1;
                    let mut end = start;
                    while end < chars.len() && !chars[end].is_whitespace() && chars[end] != ',' && chars[end] != ')' {
                        end += 1;
                    }
                    let word: String = chars[start..end].iter().collect();
                    let v = parse_date(&word).ok_or_else(|| SpiceError::TextKernel {
                        path: path.to_string(),
                        line: lineno,
                        reason: format!("cannot parse date '@{}'", word),
                    })?;
                    toks.push((lineno, Tok::Num(v)));
                    k = end;
                }
                _ => {
                    let start = k;
                    let mut end = k;
                    while end < chars.len() {
                        let ch = chars[end];
                        if ch.is_whitespace() || ch == ',' || ch == '(' || ch == ')' || ch == '=' || ch == '\'' {
                            break;
                        }
                        if ch == '+' && end + 1 < chars.len() && chars[end + 1] == '=' && end > start {
                            break;
                        }
                        end += 1;
                    }
                    let word: String = chars[start..end].iter().collect();
                    k = end;
                    // Numbers vs names: a name is followed (eventually) by '=' / '+='.
                    if let Some(v) = parse_number(&word) {
                        // Could still be a variable name that looks numeric (rare); decide by lookahead.
                        let mut j = k;
                        while j < chars.len() && chars[j].is_whitespace() {
                            j += 1;
                        }
                        let followed_by_assign = j < chars.len()
                            && (chars[j] == '=' || (chars[j] == '+' && j + 1 < chars.len() && chars[j + 1] == '='));
                        if followed_by_assign {
                            toks.push((lineno, Tok::Name(word)));
                        } else {
                            toks.push((lineno, Tok::Num(v)));
                        }
                    } else {
                        toks.push((lineno, Tok::Name(word)));
                    }
                }
            }
        }
    }
    Ok(toks)
}

fn parse_number(s: &str) -> Option<f64> {
    let t = s.replace(['D', 'd'], "E");
    t.parse::<f64>().ok()
}

/// Parse an `@`-date into seconds past J2000 (formal calendar, no leap seconds), as SPICE does
/// for kernel pool date values.
fn parse_date(s: &str) -> Option<f64> {
    const MONTHS: [&str; 12] = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"];
    let up = s.to_uppercase();
    let tpos = up
        .char_indices()
        .find(|&(i, c)| c == 'T' && i >= 8 && up[..i].ends_with(|d: char| d.is_ascii_digit()))
        .map(|(i, _)| i);
    let (date, time) = if let Some(p) = up.find('/') {
        (up[..p].to_string(), up[p + 1..].to_string())
    } else if let Some(p) = tpos {
        (up[..p].to_string(), up[p + 1..].to_string())
    } else {
        (up.clone(), String::new())
    };
    let parts: Vec<&str> = date.split('-').collect();
    if parts.len() != 3 {
        return None;
    }
    let year: i64 = parts[0].parse().ok().filter(|y: &i64| (-1_000_000..=1_000_000).contains(y))?;
    let month: i64 = match parts[1].parse::<i64>() {
        Ok(m) => m,
        Err(_) => MONTHS.iter().position(|m| parts[1].starts_with(m))? as i64 + 1,
    };
    let day: i64 = parts[2].parse().ok().filter(|d: &i64| (0..=366).contains(d))?;
    if !(1..=12).contains(&month) {
        return None;
    }
    let mut secs = 0.0;
    if !time.is_empty() {
        let tp: Vec<&str> = time.split(':').collect();
        let h: f64 = tp.first()?.parse().ok()?;
        let m: f64 = tp.get(1).map(|x| x.parse().ok()).unwrap_or(Some(0.0))?;
        let sec: f64 = tp.get(2).map(|x| x.parse().ok()).unwrap_or(Some(0.0))?;
        secs = h * 3600.0 + m * 60.0 + sec;
    }
    // Julian day number (Gregorian calendar) at noon.
    let a = (14 - month) / 12;
    let y = year + 4800 - a;
    let m = month + 12 * a - 3;
    let jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    let days_from_j2000 = (jdn - 2_451_545) as f64 - 0.5; // midnight of that date
    Some(days_from_j2000 * 86_400.0 + secs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_assignments() {
        let text = "junk\n\\begindata\nA = 1.5D3\nB = ( 1, 2\n 3 )\nC = 'x''y'\nB += 4\nNAIF_BODY_NAME += ( 'FOO' 'BAR' )\n\\begintext\nD = 5\n";
        let mut p = KernelPool::new();
        p.load_str(text, "t").unwrap();
        assert_eq!(p.get_f64("A"), Some(1500.0));
        assert_eq!(p.get_f64s("B").unwrap(), &[1.0, 2.0, 3.0, 4.0]);
        assert_eq!(p.get_str("C"), Some("x'y"));
        assert_eq!(p.get_strs("NAIF_BODY_NAME").unwrap().len(), 2);
        assert!(!p.contains("D"));
    }

    #[test]
    fn continuation_strings() {
        let text = "\\begindata\nK = ( '/a/b+' 'c.bsp' 'd.tls' )\n";
        let mut p = KernelPool::new();
        p.load_str(text, "t").unwrap();
        assert_eq!(p.get_strs_joined("K").unwrap(), vec!["/a/bc.bsp".to_string(), "d.tls".to_string()]);
    }

    #[test]
    fn markers_count_only_alone_on_their_line() {
        // pck00010.tpc's comments mention the markers in running text, with apostrophes after.
        let text = "KPL/PCK\n     \\begindata token. In order to be recognized, it doesn't matter\n\\begindata\nA = 1\n  \\begintext  \nB = 'x\n";
        let mut p = KernelPool::new();
        p.load_str(text, "t").unwrap();
        assert_eq!(p.get_f64s("A"), Some(&[1.0][..]));
        assert!(!p.contains("B"));
    }

    #[test]
    fn dates() {
        assert_eq!(parse_date("2000-JAN-01/12:00:00"), Some(0.0));
        assert_eq!(parse_date("2000-01-01T12:00:00"), Some(0.0));
        assert_eq!(parse_date("1972-JAN-1"), Some(-883_656_000.0));
    }
}
