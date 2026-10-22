use std::sync::{Arc, Mutex};
use std::time::Instant;

use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::gpu::context::{GpuContext, GpuError, Params, Precision};
use crate::gpu::wisdom_holman::GpuStats;
use crate::nbody::forces::Force;
use crate::nbody::integrators::ias15::{C, D, H, RR, SAFETY_FACTOR};
use crate::nbody::integrators::{Integrator, IAS15};
use crate::state::State;

/// Reals per test particle in the GPU state (see `ias15.wgsl`).
const ST: usize = 51;

/// IAS15 with the test particles on the GPU, each with its own adaptive timestep.
///
/// Each [`Integrator::step`] covers a fixed interval, [`GpuIAS15::timestep`] (a block: say a
/// year or a century). The massive bodies are integrated through it on the CPU with
/// [`IAS15`], keeping every step's Gauss–Radau coefficients, which give their positions anywhere
/// in the step to the step's own accuracy. Then each test particle runs IAS15 through the
/// block in its own GPU thread, in the field of the massive bodies at those positions, landing
/// on the block's end.
///
/// The test particles' steps are chosen from their own coefficients only, so a particle having
/// a close encounter doesn't shorten anyone else's steps (on the CPU, IAS15's one timestep
/// is the shortest any particle needs). Results agree with the CPU [`IAS15`] of the whole
/// system to the integrators' accuracy, not to round-off. With [`Precision::DF64`] the
/// round-off is about 2^-48 per operation; IAS15 is a 15th-order method whose error target
/// (`epsilon`, 1e-9 by default) is far above that, but plain f32 can't hold it.
///
/// Between calls the test particles stay on the GPU (with their timesteps and predicted
/// coefficients) as long as the simulation hands back the states this integrator returned.
/// Calls with forces other than Newtonian gravity, with fewer than
/// [`GpuIAS15::min_tests`] test particles or with a negative timestep run on the CPU.
pub struct GpuIAS15 {
    /// The interval each step covers, in days.
    pub timestep: f64,
    /// IAS15's error target, for the massive bodies and the test particles.
    pub epsilon: f64,
    pub precision: Precision,
    /// Fewest test particles for the GPU to be used.
    pub min_tests: usize,
    ctx: Arc<GpuContext>,
    kernels: Arc<Mutex<Vec<Kernel>>>,
    /// The massive bodies' integrator, and their states when it last stopped.
    planets: IAS15,
    planets_end: Option<Vec<State>>,
    /// For calls that run on the CPU.
    cpu: IAS15,
    resident: Option<Resident>,
    stats: Arc<Mutex<GpuStats>>,
}

impl Clone for GpuIAS15 {
    /// A copy that shares the GPU, the compiled kernels and the stats, but not the particles
    /// kept on the GPU.
    fn clone(&self) -> Self {
        GpuIAS15 {
            timestep: self.timestep,
            epsilon: self.epsilon,
            precision: self.precision,
            min_tests: self.min_tests,
            ctx: self.ctx.clone(),
            kernels: self.kernels.clone(),
            planets: self.planets.clone(),
            planets_end: self.planets_end.clone(),
            cpu: self.cpu.clone(),
            resident: None,
            stats: self.stats.clone(),
        }
    }
}

impl std::fmt::Debug for GpuIAS15 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GpuIAS15({}, {:?}, {:?})", self.timestep, self.precision, self.ctx)
    }
}

struct Kernel {
    precision: Precision,
    epsilon: f64,
    run: wgpu::ComputePipeline,
    begin: wgpu::ComputePipeline,
}

/// Test particles kept on the GPU between calls.
struct Resident {
    st: wgpu::Buffer,
    status: wgpu::Buffer,
    /// The test particles' states as the last call returned them.
    out: Vec<State>,
}

impl GpuIAS15 {
    /// A GPU IAS15 taking steps of `timestep` days. Fails if `precision` isn't available on
    /// `ctx` or its arithmetic fails the self-test there.
    pub fn new(timestep: f64, ctx: Arc<GpuContext>, precision: Precision) -> Result<GpuIAS15, GpuError> {
        ctx.check(precision)?;
        let ias = GpuIAS15 {
            timestep,
            epsilon: 1e-9,
            precision,
            min_tests: 1,
            ctx,
            kernels: Arc::new(Mutex::new(Vec::new())),
            planets: IAS15::new(timestep.min(10.0)),
            planets_end: None,
            cpu: IAS15::new(timestep.min(10.0)),
            resident: None,
            stats: Arc::new(Mutex::new(GpuStats::default())),
        };
        // Compile the kernels now rather than in the first step.
        ias.kernels();
        Ok(ias)
    }

    /// Where the time went, summed over the calls of this integrator and its clones.
    pub fn stats(&self) -> GpuStats {
        *self.stats.lock().unwrap()
    }

    pub fn reset_stats(&self) {
        *self.stats.lock().unwrap() = GpuStats::default();
    }

    fn kernels(&self) -> (wgpu::ComputePipeline, wgpu::ComputePipeline) {
        let mut kernels = self.kernels.lock().unwrap();
        if let Some(k) = kernels.iter().find(|k| k.precision == self.precision && k.epsilon == self.epsilon) {
            return (k.run.clone(), k.begin.clone());
        }
        let p = self.precision;
        let mut c = String::new();
        c.push_str(&p.constant_array("H", &H));
        c.push_str(&p.constant_array("RRI", &RR.map(|r| 1.0 / r)));
        c.push_str(&p.constant_array("C", &C));
        c.push_str(&p.constant_array("D", &D));
        // Factors of the position and velocity series at each substep (see `IAS15::step`).
        let mut xf = vec![];
        let mut vf = vec![];
        for &h in &H[1..] {
            xf.extend([7.0 * h / 9.0, 3.0 * h / 4.0, 5.0 * h / 7.0, 2.0 * h / 3.0, 3.0 * h / 5.0, h / 2.0, h / 3.0]);
            vf.extend([7.0 * h / 8.0, 6.0 * h / 7.0, 5.0 * h / 6.0, 4.0 * h / 5.0, 3.0 * h / 4.0, 2.0 * h / 3.0, h / 2.0]);
        }
        c.push_str(&p.constant_array("XF", &xf));
        c.push_str(&p.constant_array("VF", &vf));
        for (name, v) in [("THIRD", 1.0 / 3.0), ("THREE_FIFTHS", 0.6), ("TWO_THIRDS", 2.0 / 3.0), ("FIVE_SEVENTHS", 5.0 / 7.0), ("SEVEN_NINTHS", 7.0 / 9.0)] {
            c.push_str(&p.constant(name, v));
        }
        for k in [3, 4, 5, 6, 7, 8, 12, 20, 30, 42, 56, 72] {
            c.push_str(&p.constant(&format!("INV{k}"), 1.0 / k as f64));
        }
        c.push_str(&format!(
            "const SAFETY: f32 = {:?};\nconst EPS_FAC: f32 = {:?};\nconst PC_TOL: f32 = {:?};\n",
            SAFETY_FACTOR as f32,
            (self.epsilon * 5040.0).powf(1.0 / 7.0) as f32,
            p.epsilon() as f32
        ));
        let src = self.ctx.source(p, &c, &[include_str!("shaders/ias15.wgsl"), include_str!("shaders/ias15_begin.wgsl")]);
        let module = self.ctx.module("ias15", &src);
        let run = self.ctx.pipeline(&module, "main");
        let begin = self.ctx.pipeline(&module, "begin");
        kernels.push(Kernel { precision: p, epsilon: self.epsilon, run: run.clone(), begin: begin.clone() });
        (run, begin)
    }

    /// Integrate everything on the CPU through `h` days.
    fn cpu_interval(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>], h: f64) {
        self.cpu.epsilon = self.epsilon;
        ias15_interval(&mut self.cpu, states, masses, t, forces, h, |_| {});
    }

    fn step_gpu(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>]) {
        let start = Instant::now();
        let ctx = self.ctx.clone();
        let p = self.precision;
        let h = self.timestep;
        let (run_k, begin_k) = self.kernels();
        let massive: Vec<usize> = (0..states.len()).filter(|&i| masses[i] != 0.0).collect();
        let tests: Vec<usize> = (0..states.len()).filter(|&i| masses[i] == 0.0).collect();
        let n = tests.len();
        let nb = massive.len();

        // The massive bodies through the block, keeping each step.
        let tic = Instant::now();
        let mut mstates: Vec<State> = massive.iter().map(|&i| states[i]).collect();
        let m_start = mstates.clone();
        let mmasses: Vec<f64> = massive.iter().map(|&i| masses[i]).collect();
        if self.planets_end.as_ref() != Some(&mstates) {
            self.planets = IAS15::new(h.min(10.0));
        }
        self.planets.epsilon = self.epsilon;
        let mut table = Vec::new();
        let mut t_rel = 0.0;
        ias15_interval(&mut self.planets, &mut mstates, &mmasses, &mut t_rel, forces, h, |ias| {
            let (t0, dt, x0, v0, a0, b) = ias.last_step().unwrap();
            table.extend([t0, 1.0 / dt]);
            for j in 0..x0.len() {
                for v in [&x0[j], &v0[j], &a0[j], &b[j].p0, &b[j].p1, &b[j].p2, &b[j].p3, &b[j].p4, &b[j].p5, &b[j].p6] {
                    table.extend_from_slice(v.as_slice());
                }
            }
        });
        self.planets_end = Some(mstates.clone());
        let nrows = table.len() / (2 + 30 * nb);
        let cpu_seconds = tic.elapsed().as_secs_f64();

        // The test particles: still on the GPU from the last call, or uploaded.
        let current: Vec<State> = tests.iter().map(|&i| states[i]).collect();
        let reuse = matches!(&self.resident, Some(r) if r.out == current);
        if !reuse {
            let dt0 = self.planets.timestep.min(h);
            let mut st = vec![0.0; ST * n];
            for (k, s) in current.iter().enumerate() {
                st[ST * k..ST * k + 6].copy_from_slice(s);
                st[ST * k + 7] = dt0;
            }
            let mut bytes = vec![];
            p.pack(&st, &mut bytes);
            self.resident = Some(Resident { st: ctx.storage(&bytes, true), status: ctx.storage_zeroed(8 * n as u64), out: vec![] });
        }
        let res = self.resident.as_ref().unwrap();
        let mut bytes = vec![];
        p.pack(&table, &mut bytes);
        let tab_buf = ctx.storage(&bytes, false);
        let mut coef = vec![h];
        coef.extend(mmasses.iter().map(|m| GRAVITATIONAL_CONSTANT * m));
        bytes.clear();
        p.pack(&coef, &mut bytes);
        let coef_buf = ctx.storage(&bytes, false);
        let running = ctx.storage_zeroed(16);
        let mut params = Params::new(n as u32);
        (params.a, params.b, params.c) = (nrows as u32, nb as u32, 1);
        let params_buf = ctx.uniform(&params);
        let buffers = [&params_buf, &res.st, &tab_buf, &coef_buf, &res.status, &running];
        let groups = (n as u32).div_ceil(64);
        ctx.dispatch(&begin_k, &ctx.bind_group_at(&begin_k, &[(0, &params_buf), (1, &res.st), (4, &res.status)]), groups);

        // Dispatch until every particle has arrived, each dispatch near 50 ms.
        let group = ctx.bind_group(&run_k, &buffers);
        let mut max_steps = 8u32;
        let tic = Instant::now();
        let mut rounds = 0;
        let mut particle_steps = 0.0;
        loop {
            params.c = max_steps;
            ctx.write(&params_buf, bytemuck::bytes_of(&params));
            ctx.write(&running, &[0u8; 16]);
            let t_dispatch = Instant::now();
            ctx.dispatch(&run_k, &group, groups);
            let counts = ctx.read(&running);
            let left = u32::from_le_bytes(counts[..4].try_into().unwrap());
            particle_steps += u32::from_le_bytes(counts[4..8].try_into().unwrap()) as f64;
            let took = t_dispatch.elapsed().as_secs_f64();
            rounds += 1;
            if left == 0 {
                break;
            }
            max_steps = ((max_steps as f64 * 0.05 / took.max(1e-4)) as u32).clamp(1, 1 << 16);
            if rounds > 1_000_000 {
                panic!("GpuIAS15: test particles not arriving");
            }
        }
        let gpu_seconds = tic.elapsed().as_secs_f64();

        let st = p.unpack(&ctx.read(&res.st));
        let status: Vec<u32> = bytemuck::cast_slice(&ctx.read(&res.status)).to_vec();
        for (k, &i) in massive.iter().enumerate() {
            states[i] = mstates[k];
        }
        let mut failed = vec![];
        for (k, &i) in tests.iter().enumerate() {
            if status[2 * k] != 1 {
                failed.push(i);
                continue;
            }
            states[i] = st[ST * k..ST * k + 6].try_into().unwrap();
        }
        if !failed.is_empty() {
            // Those particles again on the CPU, with the massive bodies.
            let keep: Vec<usize> = massive.iter().copied().chain(failed.iter().copied()).collect();
            let mut sub: Vec<State> = m_start.clone();
            sub.extend(failed.iter().map(|&i| current[tests.iter().position(|&j| j == i).unwrap()]));
            let sub_m: Vec<f64> = keep.iter().map(|&i| masses[i]).collect();
            let mut t_sub = 0.0;
            self.cpu = IAS15::new(h.min(10.0));
            self.cpu_interval(&mut sub, &sub_m, &mut t_sub, forces, h);
            for (s, &i) in sub.iter().zip(&keep) {
                if masses[i] == 0.0 {
                    states[i] = *s;
                }
            }
            // The GPU copy is stale now.
            self.resident = None;
        }
        if let Some(r) = self.resident.as_mut() {
            r.out = tests.iter().map(|&i| states[i]).collect();
        }
        *t += h;

        let mut stats = self.stats.lock().unwrap();
        stats.particle_steps += particle_steps;
        stats.cpu_seconds += cpu_seconds;
        stats.gpu_wait_seconds += gpu_seconds;
        stats.total_seconds += start.elapsed().as_secs_f64();
        stats.cpu_fallbacks += failed.len();
    }
}

/// Take IAS15 steps through `h` days, landing on the end (as `Simulation::integrate` does), and
/// call `each` after every step.
fn ias15_interval(ias: &mut IAS15, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>], h: f64, mut each: impl FnMut(&IAS15)) {
    let end = *t + h;
    loop {
        let rem = end - *t;
        if rem.abs() <= 1e-14 * h.abs().max(1.0) {
            break;
        }
        let full = ias.timestep();
        if full.abs() > rem.abs() || full.signum() != rem.signum() {
            let start = *t;
            ias.set_timestep(rem);
            ias.step(states, masses, t, forces);
            each(ias);
            if ((*t - start) - rem).abs() < 0.5 * rem.abs() {
                ias.set_timestep(full.abs() * rem.signum());
            }
        } else {
            ias.step(states, masses, t, forces);
            each(ias);
        }
    }
}

impl Integrator for GpuIAS15 {
    fn step(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>]) {
        let newtonian_only = forces.len() == 1 && forces[0].is_newtonian_gravity();
        let n_tests = masses.iter().filter(|&&m| m == 0.0).count();
        let any_massive = masses.iter().any(|&m| m != 0.0);
        if newtonian_only && any_massive && n_tests >= self.min_tests.max(1) && self.timestep > 0.0 {
            self.step_gpu(states, masses, t, forces);
        } else {
            self.resident = None;
            let h = self.timestep;
            self.cpu_interval(states, masses, t, forces, h);
        }
    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }

    fn fixed_timestep(&self) -> bool {
        true
    }
}
