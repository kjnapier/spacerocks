//! The small localhost web server behind spacerocks' live pages: the simulation monitor
//! ([`crate::nbody::Simulation::monitor`]) and the orbit-fit dashboard
//! ([`crate::orbfit::FitDashboard`]).
//!
//! Standard library only. Each page binds `127.0.0.1`, answers one request per connection on
//! its own thread, and stops once the state it serves is dropped. Besides the page's own
//! routes, every server answers the shared retro assets: `font.ttf` (a 5x7 pixel font built by
//! `scripts/monitor_font.py`), `retro.css` (the 8-bit look: palette, RPG windows, buttons and
//! toggles) and `retro.js` (formatting helpers and pixel drawing on low-resolution canvases).
//!
//! Pages change things with `POST` requests, which are only accepted with an `X-Spacerocks`
//! header and a local `Host`, so other web pages open in the browser can't send them.
use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::Weak;
use std::time::Duration;

/// The pixel font, built by `scripts/monitor_font.py`.
const FONT: &[u8] = include_bytes!("font.ttf");
const RETRO_CSS: &str = include_str!("retro.css");
const RETRO_JS: &str = include_str!("retro.js");

/// A request's method, path and query string.
pub(crate) struct Request {
    pub method: String,
    pub path: String,
    pub query: String,
}

impl Request {
    /// The decoded value of `key` in the query string.
    pub fn param(&self, key: &str) -> Option<String> {
        self.query.split('&').filter_map(|kv| kv.split_once('=')).find(|(k, _)| *k == key).map(|(_, v)| percent_decode(v))
    }
}

/// A response: status line (like `200 OK`), content type and body.
pub(crate) struct Response {
    pub status: &'static str,
    pub content_type: &'static str,
    pub body: Vec<u8>,
}

impl Response {
    pub fn ok(content_type: &'static str, body: impl Into<Vec<u8>>) -> Response {
        Response { status: "200 OK", content_type, body: body.into() }
    }

    pub fn json<T: serde::Serialize>(value: &T) -> Response {
        Response::ok("application/json", serde_json::to_vec(value).unwrap_or_default())
    }

    pub fn html(page: &str) -> Response {
        Response::ok("text/html; charset=utf-8", page)
    }

    pub fn text(status: &'static str, body: impl Into<String>) -> Response {
        Response { status, content_type: "text/plain", body: body.into().into_bytes() }
    }

    pub fn not_found() -> Response {
        Response::text("404 Not Found", "not found")
    }

    /// `ok` for `Ok`, 400 with the message for `Err`.
    pub fn outcome(r: Result<(), String>) -> Response {
        match r {
            Ok(()) => Response::text("200 OK", "ok"),
            Err(e) => Response::text("400 Bad Request", e),
        }
    }
}

/// Bind `127.0.0.1:port` (0 picks a free port) and answer requests with `handle` on a thread
/// named `name` until `state` is dropped. Returns the page's URL.
pub(crate) fn start<S, F>(port: u16, name: &str, state: Weak<S>, handle: F) -> std::io::Result<String>
where
    S: Send + Sync + 'static,
    F: Fn(&Request, &S) -> Response + Send + 'static,
{
    let listener = TcpListener::bind(("127.0.0.1", port))?;
    listener.set_nonblocking(true)?;
    let port = listener.local_addr()?.port();
    std::thread::Builder::new().name(name.into()).spawn(move || serve(listener, state, port, handle))?;
    Ok(format!("http://127.0.0.1:{port}/"))
}

fn serve<S, F>(listener: TcpListener, state: Weak<S>, port: u16, handle: F)
where
    F: Fn(&Request, &S) -> Response,
{
    loop {
        if state.strong_count() == 0 {
            return;
        }
        match listener.accept() {
            Ok((stream, _)) => {
                let Some(state) = state.upgrade() else { return };
                let _ = respond(stream, port, |r| handle(r, &state));
            }
            Err(_) => std::thread::sleep(Duration::from_millis(25)),
        }
    }
}

fn respond(mut stream: TcpStream, port: u16, handle: impl Fn(&Request) -> Response) -> std::io::Result<()> {
    stream.set_nonblocking(false)?;
    stream.set_read_timeout(Some(Duration::from_secs(2)))?;
    // Read up to the end of the headers (requests carry no body), which may come in several
    // pieces: closing the connection with part of a request unread would reset it.
    let mut buf = Vec::with_capacity(1024);
    let mut chunk = [0u8; 1024];
    while !buf.windows(4).any(|w| w == b"\r\n\r\n") && buf.len() < 16 * 1024 {
        let n = stream.read(&mut chunk)?;
        if n == 0 {
            break;
        }
        buf.extend_from_slice(&chunk[..n]);
    }
    let head = String::from_utf8_lossy(&buf).into_owned();
    let mut words = head.split_whitespace();
    let method = words.next().unwrap_or("GET").to_string();
    let target = words.next().unwrap_or("/");
    let (path, query) = target.split_once('?').unwrap_or((target, ""));
    let request = Request { method, path: path.to_string(), query: query.to_string() };

    let response = match (request.method.as_str(), request.path.as_str()) {
        ("GET", "/font.ttf") => Response::ok("font/ttf", FONT),
        ("GET", "/retro.css") => Response::ok("text/css; charset=utf-8", RETRO_CSS),
        ("GET", "/retro.js") => Response::ok("text/javascript; charset=utf-8", RETRO_JS),
        ("POST", _) if !from_page(&head, port) => Response::text("403 Forbidden", "forbidden"),
        _ => handle(&request),
    };
    write!(
        stream,
        "HTTP/1.1 {}\r\nContent-Type: {}\r\nContent-Length: {}\r\nCache-Control: no-store\r\nConnection: close\r\n\r\n",
        response.status,
        response.content_type,
        response.body.len()
    )?;
    stream.write_all(&response.body)?;
    stream.flush()
}

/// Whether a request came from the page itself: it carries the `X-Spacerocks` header (which a
/// page from another site can't add without a CORS preflight these servers never answer) and a
/// local `Host` (so a site rebinding its name to 127.0.0.1 is turned away).
fn from_page(head: &str, port: u16) -> bool {
    let mut marked = false;
    let mut local = false;
    for line in head.lines().skip(1) {
        let Some((name, value)) = line.split_once(':') else { continue };
        let value = value.trim();
        if name.eq_ignore_ascii_case("x-spacerocks") {
            marked = true;
        } else if name.eq_ignore_ascii_case("host") {
            local = value == format!("127.0.0.1:{port}") || value == format!("localhost:{port}");
        }
    }
    marked && local
}

/// Decode a query-string value (`+` and `%XX`).
pub(crate) fn percent_decode(v: &str) -> String {
    let hex = |b: u8| (b as char).to_digit(16).map(|d| d as u8);
    let bytes = v.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        match (bytes[i], bytes.get(i + 1).copied().and_then(hex), bytes.get(i + 2).copied().and_then(hex)) {
            (b'%', Some(hi), Some(lo)) => {
                out.push(hi * 16 + lo);
                i += 3;
            }
            (b'+', ..) => {
                out.push(b' ');
                i += 1;
            }
            (b, ..) => {
                out.push(b);
                i += 1;
            }
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}
