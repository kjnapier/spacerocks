use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use wgpu::util::DeviceExt;

/// How the GPU kernels represent real numbers.
///
/// WGSL has no f64 in its core: native f64 needs the adapter's `SHADER_F64` feature (Vulkan and
/// DX12 drivers that support it; never Metal, so never Apple GPUs), and consumer GPUs run it at
/// 1/16 to 1/64 of their f32 rate. `DF64` emulates about 48 bits of mantissa with pairs of f32s
/// and works on every adapter. Plain `F32` is here to measure how wrong it goes; it isn't
/// accurate enough to integrate with (see `notes` in the module docs).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Precision {
    F32,
    DF64,
    F64,
}

impl Precision {
    /// Bytes per real number in a buffer.
    pub fn bytes(self) -> usize {
        match self {
            Precision::F32 => 4,
            Precision::DF64 | Precision::F64 => 8,
        }
    }

    /// Unit round-off: the relative error of one rounding.
    pub fn epsilon(self) -> f64 {
        match self {
            Precision::F32 => 2f64.powi(-24),
            Precision::DF64 => 2f64.powi(-48),
            Precision::F64 => 2f64.powi(-53),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Precision::F32 => "f32",
            Precision::DF64 => "df64",
            Precision::F64 => "f64",
        }
    }

    pub fn parse(s: &str) -> Option<Precision> {
        match s.to_ascii_lowercase().as_str() {
            "f32" => Some(Precision::F32),
            "df64" => Some(Precision::DF64),
            "f64" => Some(Precision::F64),
            _ => None,
        }
    }

    /// Append `values` to `out` in this precision's buffer layout.
    pub fn pack(self, values: &[f64], out: &mut Vec<u8>) {
        for &v in values {
            match self {
                Precision::F32 => out.extend_from_slice(&(v as f32).to_le_bytes()),
                Precision::DF64 => {
                    let [hi, lo] = split(v);
                    out.extend_from_slice(&hi.to_le_bytes());
                    out.extend_from_slice(&lo.to_le_bytes());
                }
                Precision::F64 => out.extend_from_slice(&v.to_le_bytes()),
            }
        }
    }

    /// Read back real numbers packed by [`Precision::pack`].
    pub fn unpack(self, bytes: &[u8]) -> Vec<f64> {
        let f = |b: &[u8]| f32::from_le_bytes(b.try_into().unwrap()) as f64;
        match self {
            Precision::F32 => bytes.chunks_exact(4).map(f).collect(),
            Precision::DF64 => bytes.chunks_exact(8).map(|b| f(&b[..4]) + f(&b[4..])).collect(),
            Precision::F64 => bytes.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect(),
        }
    }

    /// The WGSL library for this precision: the type `R` and its arithmetic.
    fn library(self) -> &'static str {
        match self {
            Precision::F32 => include_str!("shaders/real_f32.wgsl"),
            Precision::DF64 => include_str!("shaders/real_df64.wgsl"),
            Precision::F64 => include_str!("shaders/real_f64.wgsl"),
        }
    }

    /// A WGSL expression of type `R` with value `v`, as close as this precision gets.
    pub(crate) fn literal(self, v: f64) -> String {
        match self {
            Precision::F32 => format!("{:?}", v as f32),
            Precision::DF64 => {
                let [hi, lo] = split(v);
                format!("vec2<f32>({hi:?}, {lo:?})")
            }
            Precision::F64 => format!("{v:e}lf"),
        }
    }

    /// A WGSL constant of type `R` with value `v`.
    pub(crate) fn constant(self, name: &str, v: f64) -> String {
        format!("const {name}: R = {};\n", self.literal(v))
    }

    /// A WGSL constant array of `R`.
    pub(crate) fn constant_array(self, name: &str, v: &[f64]) -> String {
        let items: Vec<String> = v.iter().map(|&x| self.literal(x)).collect();
        format!("const {name} = array<R, {}>({});\n", v.len(), items.join(", "))
    }
}

/// `v` as an f32 pair `hi + lo`.
fn split(v: f64) -> [f32; 2] {
    let hi = v as f32;
    let lo = (v - hi as f64) as f32;
    [hi, lo]
}

#[derive(Debug, thiserror::Error)]
pub enum GpuError {
    #[error("no GPU adapter: {0}")]
    NoAdapter(String),
    #[error("could not open the GPU: {0}")]
    Device(String),
    #[error("{0} is not supported on {1}")]
    Unsupported(&'static str, String),
    #[error("{0} arithmetic fails its self-test on {1}: {2}")]
    SelfTest(&'static str, String, String),
    #[error("GPU error: {0}")]
    Other(String),
}

/// An open GPU: the device the integrators run their kernels on.
///
/// The adapter is the system's high-performance one, unless `SPACEROCKS_GPU` names (part of) the
/// name of another: `SPACEROCKS_GPU=llvmpipe` picks Mesa's software Vulkan driver, which is what
/// the tests run on when there is no GPU. `WGPU_BACKEND` limits the backends (`vulkan`, `metal`,
/// `dx12`, `gl`).
pub struct GpuContext {
    pub(crate) device: wgpu::Device,
    pub(crate) queue: wgpu::Queue,
    info: wgpu::AdapterInfo,
    f64: bool,
    /// Results of the arithmetic self-test per precision: the largest relative error.
    checked: Mutex<HashMap<Precision, Result<f64, String>>>,
    /// Whether fma is fused (exact) in shaders, once the self-test has found out.
    fma_exact: Mutex<Option<bool>>,
}

impl std::fmt::Debug for GpuContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GpuContext({})", self.describe())
    }
}

impl GpuContext {
    /// Open the GPU.
    pub fn new() -> Result<Arc<GpuContext>, GpuError> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle_from_env());
        let wanted = std::env::var("SPACEROCKS_GPU").ok().map(|s| s.to_lowercase());
        let adapter = match &wanted {
            Some(name) => {
                let adapters = pollster::block_on(instance.enumerate_adapters(wgpu::Backends::all()));
                let names: Vec<String> = adapters.iter().map(|a| a.get_info().name).collect();
                adapters.into_iter().find(|a| a.get_info().name.to_lowercase().contains(name.as_str()))
                    .ok_or_else(|| GpuError::NoAdapter(format!("none matches SPACEROCKS_GPU={name} (found {names:?})")))?
            }
            None => pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                ..Default::default()
            }))
            .map_err(|e| GpuError::NoAdapter(e.to_string()))?,
        };
        let info = adapter.get_info();
        let f64 = adapter.features().contains(wgpu::Features::SHADER_F64);
        let features = if f64 { wgpu::Features::SHADER_F64 } else { wgpu::Features::empty() };
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("spacerocks"),
            required_features: features,
            required_limits: adapter.limits(),
            ..Default::default()
        }))
        .map_err(|e| GpuError::Device(e.to_string()))?;
        Ok(Arc::new(GpuContext { device, queue, info, f64, checked: Mutex::new(HashMap::new()), fma_exact: Mutex::new(None) }))
    }

    /// The adapter's name and backend.
    pub fn describe(&self) -> String {
        format!("{} ({:?}, {:?})", self.info.name, self.info.backend, self.info.device_type)
    }

    /// Whether fma is exact in this adapter's shaders (known after checking `DF64`).
    pub fn fma_exact(&self) -> Option<bool> {
        *self.fma_exact.lock().unwrap()
    }

    /// Whether the adapter has native f64 in shaders.
    pub fn has_f64(&self) -> bool {
        self.f64
    }

    /// Check that `precision` works here: available, and its arithmetic right (for `DF64`, that
    /// the shader compiler kept the error-free transformations exact). Runs a self-test the
    /// first time and returns the largest relative error it saw.
    pub fn check(&self, precision: Precision) -> Result<f64, GpuError> {
        if precision == Precision::F64 && !self.f64 {
            return Err(GpuError::Unsupported("f64", self.describe()));
        }
        let mut checked = self.checked.lock().unwrap();
        if precision == Precision::DF64 && !checked.contains_key(&precision) {
            // Exact products with fma if the driver fuses it, with Dekker's splitting if not.
            *self.fma_exact.lock().unwrap() = Some(true);
            if self.self_test(precision).is_err() {
                *self.fma_exact.lock().unwrap() = Some(false);
            }
        }
        let result = checked.entry(precision).or_insert_with(|| self.self_test(precision)).clone();
        result.map_err(|e| GpuError::SelfTest(precision.name(), self.describe(), e))
    }

    /// Run add, multiply, divide and square root on awkward operands and compare with f64.
    fn self_test(&self, precision: Precision) -> Result<f64, String> {
        // Operands near each other (cancellation), far apart, and of mixed signs, as f64 values
        // that the precision represents exactly, so only the operations' errors show.
        let mut seed: u64 = 7;
        let mut rnd = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut a = Vec::new();
        let mut b = Vec::new();
        for k in 0..4096 {
            let x = (rnd() - 0.5) * 10f64.powf(8.0 * rnd() - 4.0);
            let y = match k % 4 {
                0 => -x * (1.0 + 1e-9 * (rnd() - 0.5)),
                1 => x * (1.0 + 1e-3 * rnd()),
                _ => (rnd() - 0.5) * 10f64.powf(8.0 * rnd() - 4.0),
            };
            let exact = |v: f64| precision.unpack(&{ let mut o = vec![]; precision.pack(&[v], &mut o); o })[0];
            a.push(exact(x));
            b.push(exact(y));
        }
        let n = a.len();
        let got = self.arithmetic(precision, &a, &b);

        // Worst relative error of each operation, and where.
        let mut worst = [0.0f64; 4];
        let mut what = vec![String::new(); 4];
        for i in 0..n {
            let want = [a[i] + b[i], a[i] * b[i], a[i] / b[i], a[i].abs().sqrt()];
            // Relative to the operands' scale for the sum (its result can cancel to nothing).
            let scale = [a[i].abs().max(b[i].abs()), want[1].abs(), want[2].abs(), want[3].abs()];
            for k in 0..4 {
                let err = (got[4 * i + k] - want[k]).abs() / scale[k];
                if !(err <= worst[k]) {
                    worst[k] = if err.is_nan() { f64::INFINITY } else { err };
                    what[k] = format!("{} of {:e} and {:e}: {:e}, want {:e}", ["sum", "product", "quotient", "square root"][k], a[i], b[i], got[4 * i + k], want[k]);
                }
            }
        }
        // Each operation should be good to a few units of round-off.
        let bad: Vec<String> = (0..4).filter(|&k| worst[k] > 16.0 * precision.epsilon()).map(|k| format!("relative error {:e} in the {}", worst[k], what[k])).collect();
        if !bad.is_empty() {
            return Err(bad.join("; "));
        }
        let worst = worst.iter().cloned().fold(0.0, f64::max);
        if precision != Precision::F32 {
            self.kepler_test(precision)?;
        }
        Ok(worst)
    }

    /// The arithmetic inside a large kernel, where the compiler has more to rearrange: Kepler
    /// steps on bound and unbound orbits against the CPU's solver.
    fn kepler_test(&self, precision: Precision) -> Result<(), String> {
        use crate::constants::GRAVITATIONAL_CONSTANT as MU;
        let mut cases = vec![];
        for k in 0..256 {
            let f = k as f64;
            let (a, e) = if k % 8 == 7 { (-2.0 - 0.01 * f, 1.5) } else { (1.0 + 0.2 * f, 0.003 * f) };
            let p = a * (1.0 - e * e);
            let nu = if e < 1.0 { 0.37 * f } else { 0.5 * (0.1 * f).sin() };
            let r = p / (1.0 + e * f64::cos(nu));
            let s = (MU / p).sqrt();
            cases.push(([r * nu.cos(), r * nu.sin(), 0.01 * r], [-s * nu.sin(), s * (e + nu.cos()), 1e-4 * s], MU, 3.0 + 7.0 * f));
        }
        let got = crate::gpu::wisdom_holman::kepler_steps(self, precision, &cases);
        for ((r0, v0, mu, dt), (g, ok)) in cases.iter().zip(&got) {
            let (r, v) = crate::transforms::universal_kepler_step(&(*r0).into(), &(*v0).into(), *mu, *dt).map_err(|e| e.to_string())?;
            let err = ((g[0] - r.x).powi(2) + (g[1] - r.y).powi(2) + (g[2] - r.z).powi(2)).sqrt() / r.norm()
                + ((g[3] - v.x).powi(2) + (g[4] - v.y).powi(2) + (g[5] - v.z).powi(2)).sqrt() / v.norm();
            if !ok || !(err < 1e5 * precision.epsilon()) {
                return Err(format!("Kepler step from {r0:?}, {v0:?} by {dt} days off by {err:e} (converged: {ok})"));
            }
        }
        Ok(())
    }

    /// `a + b`, `a * b`, `a / b` and `sqrt(|a|)` for each pair, on the GPU in `precision`.
    #[doc(hidden)]
    pub fn arithmetic(&self, precision: Precision, a: &[f64], b: &[f64]) -> Vec<f64> {
        let n = a.len();
        let src = self.source(precision, "", &[include_str!("shaders/selftest.wgsl")]);
        let module = self.module("selftest", &src);
        let pipeline = self.pipeline(&module, "main");
        let mut bytes = vec![];
        precision.pack(a, &mut bytes);
        let buf_a = self.storage(&bytes, false);
        bytes.clear();
        precision.pack(b, &mut bytes);
        let buf_b = self.storage(&bytes, false);
        let out = self.storage_zeroed((4 * n * precision.bytes()) as u64);
        let params = self.uniform(&Params::new(n as u32));
        let group = self.bind_group(&pipeline, &[&params, &buf_a, &buf_b, &out]);
        self.dispatch(&pipeline, &group, (n as u32).div_ceil(64));
        precision.unpack(&self.read(&out))
    }

    /// Shader source for `precision`: its arithmetic, the vector helpers, `extra` (constant
    /// declarations) and `parts`.
    pub(crate) fn source(&self, precision: Precision, extra: &str, parts: &[&str]) -> String {
        let mut s = String::new();
        if precision == Precision::DF64 {
            s.push_str(&format!("const FMA_EXACT: bool = {};\n", self.fma_exact.lock().unwrap().unwrap_or(false)));
        }
        s.push_str(precision.library());
        s.push_str(include_str!("shaders/vec.wgsl"));
        s.push_str(extra);
        for p in parts {
            s.push_str(p);
        }
        s
    }

    pub(crate) fn module(&self, label: &str, source: &str) -> wgpu::ShaderModule {
        self.device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some(label), source: wgpu::ShaderSource::Wgsl(source.into()) })
    }

    pub(crate) fn pipeline(&self, module: &wgpu::ShaderModule, entry: &str) -> wgpu::ComputePipeline {
        self.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(entry),
            layout: None,
            module,
            entry_point: Some(entry),
            compilation_options: Default::default(),
            cache: None,
        })
    }

    /// A storage buffer holding `bytes` (at least 16 bytes long, as wgpu wants).
    pub(crate) fn storage(&self, bytes: &[u8], writable: bool) -> wgpu::Buffer {
        let mut usage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST;
        if writable {
            usage |= wgpu::BufferUsages::COPY_SRC;
        }
        let mut padded;
        let contents = if bytes.len() < 16 {
            padded = bytes.to_vec();
            padded.resize(16, 0);
            &padded[..]
        } else {
            bytes
        };
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: None, contents, usage })
    }

    pub(crate) fn storage_zeroed(&self, size: u64) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: size.max(16),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        })
    }

    pub(crate) fn uniform<T: bytemuck::Pod>(&self, value: &T) -> wgpu::Buffer {
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: None,
            contents: bytemuck::bytes_of(value),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        })
    }

    pub(crate) fn bind_group(&self, pipeline: &wgpu::ComputePipeline, buffers: &[&wgpu::Buffer]) -> wgpu::BindGroup {
        let all: Vec<(u32, &wgpu::Buffer)> = buffers.iter().enumerate().map(|(i, b)| (i as u32, *b)).collect();
        self.bind_group_at(pipeline, &all)
    }

    /// A bind group with the given bindings only (an entry point's layout has just the ones it
    /// uses).
    pub(crate) fn bind_group_at(&self, pipeline: &wgpu::ComputePipeline, buffers: &[(u32, &wgpu::Buffer)]) -> wgpu::BindGroup {
        let entries: Vec<wgpu::BindGroupEntry> = buffers.iter()
            .map(|(i, b)| wgpu::BindGroupEntry { binding: *i, resource: b.as_entire_binding() })
            .collect();
        self.device.create_bind_group(&wgpu::BindGroupDescriptor { label: None, layout: &pipeline.get_bind_group_layout(0), entries: &entries })
    }

    /// Run `groups` workgroups of `pipeline` and wait for them.
    pub(crate) fn dispatch(&self, pipeline: &wgpu::ComputePipeline, group: &wgpu::BindGroup, groups: u32) {
        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(pipeline);
            pass.set_bind_group(0, group, &[]);
            pass.dispatch_workgroups(groups.max(1), 1, 1);
        }
        self.queue.submit([encoder.finish()]);
        self.wait();
    }

    pub(crate) fn wait(&self) {
        self.device.poll(wgpu::PollType::wait_indefinitely()).expect("GPU poll failed");
    }

    /// Copy `buffer` back to the host.
    pub(crate) fn read(&self, buffer: &wgpu::Buffer) -> Vec<u8> {
        let size = buffer.size();
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut encoder = self.device.create_command_encoder(&Default::default());
        encoder.copy_buffer_to_buffer(buffer, 0, &staging, 0, size);
        self.queue.submit([encoder.finish()]);
        let slice = staging.slice(..);
        slice.map_async(wgpu::MapMode::Read, |r| r.expect("GPU buffer map failed"));
        self.wait();
        let bytes = slice.get_mapped_range().expect("GPU buffer map failed").to_vec();
        staging.unmap();
        bytes
    }

    pub(crate) fn write(&self, buffer: &wgpu::Buffer, bytes: &[u8]) {
        self.queue.write_buffer(buffer, 0, bytes);
    }
}

/// The uniform block every kernel starts with. `mask` is all ones, read at run time so the
/// compiler can't fold it (see `real_df64.wgsl`); the rest is per kernel.
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub(crate) struct Params {
    pub n: u32,
    pub mask: u32,
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub d: u32,
    pub e: u32,
    pub f: u32,
}

impl Params {
    pub fn new(n: u32) -> Params {
        Params { n, mask: u32::MAX, a: 0, b: 0, c: 0, d: 0, e: 0, f: 0 }
    }
}
