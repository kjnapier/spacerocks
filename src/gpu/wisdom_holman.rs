use std::sync::{Arc, Mutex};
use std::time::Instant;

use nalgebra::Vector3;

use crate::gpu::context::{GpuContext, GpuError, Params, Precision};
use crate::nbody::forces::{central_body, Force};
use crate::nbody::integrators::wisdom_holman::MassiveRun;
use crate::nbody::integrators::{Integrator, WisdomHolman};
use crate::state::{from_pv, pv, State};

/// [`WisdomHolman`] with the test particles on the GPU.
///
/// [`Integrator::steps`] runs the massive bodies on the CPU, a block of steps at a time, exactly
/// as [`WisdomHolman::steps`] does, and a kernel takes every test particle through each block,
/// one GPU thread per particle. The test particles stay on the GPU in democratic heliocentric
/// coordinates for the whole call, so a call should cover many steps (the simulation's
/// `integrate` hands over all of them at once).
///
/// With [`Precision::F64`] the results agree with [`WisdomHolman`] to round-off (they aren't
/// bit-identical: the GPU skips the round trip to inertial coordinates between steps). With
/// [`Precision::DF64`] each step carries about 2^-48 relative error instead of 2^-53. Calls with
/// forces other than Newtonian gravity, without test particles or with fewer than
/// [`GpuWisdomHolman::min_work`] test-particle steps, and single [`Integrator::step`]s, run on
/// the CPU.
///
/// A test particle whose Kepler solve fails on the GPU (it never has in the tests) is redone on
/// the CPU for the whole call.
#[derive(Clone)]
pub struct GpuWisdomHolman {
    /// Fixed timestep, in days.
    pub timestep: f64,
    pub precision: Precision,
    /// Compute the kicks in f32 (positions differences and the state stay in `precision`). The
    /// kicks are small corrections to the velocities, so their relative error matters less; see
    /// the precision study in the module docs.
    pub fast_kick: bool,
    /// Steps per block of the massive bodies (each block's record is uploaded at once).
    pub block: usize,
    /// Fewest test-particle steps in a call for the GPU to be used.
    pub min_work: usize,
    cpu: WisdomHolman,
    ctx: Arc<GpuContext>,
    kernels: Arc<Mutex<Vec<Kernel>>>,
    stats: Arc<Mutex<GpuStats>>,
}

/// Where a GPU integrator's time went, summed over its calls.
#[derive(Clone, Copy, Debug, Default)]
pub struct GpuStats {
    /// Test-particle steps taken on the GPU.
    pub particle_steps: f64,
    /// Seconds stepping the massive bodies on the CPU.
    pub cpu_seconds: f64,
    /// Seconds waiting for the GPU (after the CPU's own work for the block).
    pub gpu_wait_seconds: f64,
    /// Seconds in the GPU calls, end to end.
    pub total_seconds: f64,
    /// Test particles redone on the CPU after a failed Kepler solve.
    pub cpu_fallbacks: usize,
}

struct Kernel {
    precision: Precision,
    fast_kick: bool,
    pipeline: wgpu::ComputePipeline,
}

impl std::fmt::Debug for GpuWisdomHolman {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GpuWisdomHolman({}, {:?}, {:?})", self.timestep, self.precision, self.ctx)
    }
}

/// Precision-dependent constants of the Kepler solver.
pub(crate) fn kepler_constants(p: Precision) -> String {
    let u = p.epsilon() as f32;
    let mut s = format!(
        "const KEPLER_K_TOL: f32 = {:?};\nconst KEPLER_DX_TOL: f32 = {:?};\nconst KEPLER_STALL: f32 = {:?};\n",
        0.1 * u,
        2.0 * u,
        (1e4 * u).min(1e-5)
    );
    let mut fact = 1.0f64;
    for k in 1..=15 {
        fact *= k as f64;
        if k >= 2 {
            s.push_str(&p.constant(&format!("INV_F{k}"), 1.0 / fact));
        }
    }
    s.push_str(&p.constant("TAU", std::f64::consts::TAU));
    s
}

impl GpuWisdomHolman {
    /// A GPU Wisdom–Holman integrator with timestep `timestep` (days). Fails if `precision`
    /// isn't available on `ctx` or its arithmetic fails the self-test there.
    pub fn new(timestep: f64, ctx: Arc<GpuContext>, precision: Precision) -> Result<GpuWisdomHolman, GpuError> {
        ctx.check(precision)?;
        Ok(GpuWisdomHolman {
            timestep,
            precision,
            fast_kick: false,
            block: 4096,
            min_work: 1 << 16,
            cpu: WisdomHolman::new(timestep),
            ctx,
            kernels: Arc::new(Mutex::new(Vec::new())),
            stats: Arc::new(Mutex::new(GpuStats::default())),
        })
    }

    pub fn context(&self) -> &Arc<GpuContext> {
        &self.ctx
    }

    pub fn stats(&self) -> GpuStats {
        *self.stats.lock().unwrap()
    }

    pub fn reset_stats(&self) {
        *self.stats.lock().unwrap() = GpuStats::default();
    }

    fn pipeline(&self) -> wgpu::ComputePipeline {
        let mut kernels = self.kernels.lock().unwrap();
        if let Some(k) = kernels.iter().find(|k| k.precision == self.precision && k.fast_kick == self.fast_kick) {
            return k.pipeline.clone();
        }
        let consts = format!("{}const FAST_KICK: bool = {};\n", kepler_constants(self.precision), self.fast_kick);
        let src = self.ctx.source(self.precision, &consts, &[include_str!("shaders/kepler.wgsl"), include_str!("shaders/wisdom_holman.wgsl")]);
        let module = self.ctx.module("wisdom_holman", &src);
        let pipeline = self.ctx.pipeline(&module, "main");
        kernels.push(Kernel { precision: self.precision, fast_kick: self.fast_kick, pipeline: pipeline.clone() });
        pipeline
    }

    /// `n` steps with the test particles on the GPU.
    fn steps_gpu(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, central: usize, forces: &[Box<dyn Force + Send + Sync>], n: usize) {
        let start = Instant::now();
        let ctx = self.ctx.clone();
        let p = self.precision;
        let h = self.timestep;
        let pipeline = self.pipeline();

        let initial: Vec<State> = states.to_vec();
        let mut run = MassiveRun::new(states, masses, central);
        let np = run.pullers.len();
        let stride = 6 + 6 * np;
        let test_idx: Vec<usize> = (0..states.len()).filter(|&i| masses[i] == 0.0).collect();
        let n_tests = test_idx.len();

        let mut cpu_seconds = 0.0;
        let mut wait_seconds = 0.0;
        let block = self.block.max(1);
        let record = |run: &mut MassiveRun, len: usize, t: &mut f64| {
            run.clear();
            for _ in 0..len {
                run.step(forces, h, |_, _| true);
                *t += h;
            }
        };

        // The first block, then the test particles in democratic heliocentric coordinates as
        // its first step loads them.
        let tic = Instant::now();
        let mut len = block.min(n);
        record(&mut run, len, t);
        cpu_seconds += tic.elapsed().as_secs_f64();
        let mut st = Vec::with_capacity(6 * n_tests);
        for &i in &test_idx {
            let (x, v) = pv(&states[i]);
            let (q, u) = run.test_load(0, &x, &v);
            st.extend_from_slice(&[q.x, q.y, q.z, u.x, u.y, u.z]);
        }
        let mut bytes = Vec::new();
        p.pack(&st, &mut bytes);
        let st_buf = ctx.storage(&bytes, true);
        let flags_buf = ctx.storage_zeroed(4 * n_tests as u64);
        let mut coef = vec![run.gm, h, 0.5 * h];
        coef.extend(run.pullers.iter().map(|m| crate::constants::GRAVITATIONAL_CONSTANT * m));
        bytes.clear();
        p.pack(&coef, &mut bytes);
        let coef_buf = ctx.storage(&bytes, false);
        let tab_buf = ctx.storage_zeroed((block * stride * p.bytes()) as u64);
        let params_buf = ctx.uniform(&Params::new(n_tests as u32));
        let group = ctx.bind_group(&pipeline, &[&params_buf, &st_buf, &tab_buf, &coef_buf, &flags_buf]);
        let groups = (n_tests as u32).div_ceil(64);

        // Steps per dispatch, adjusted to keep each near 50 ms (GPU watchdogs stop long ones).
        let mut per_dispatch = 64usize;
        let mut done = 0;
        let mut last_frame = run.frames[0];
        let mut table = Vec::with_capacity(block * stride);
        loop {
            // Upload this block's record and queue its dispatches.
            table.clear();
            for (k, f) in run.frames.iter().enumerate() {
                table.extend_from_slice(f.dq1.as_slice());
                table.extend_from_slice(f.dq2.as_slice());
                for q in &run.kicks[2 * np * k..2 * np * (k + 1)] {
                    table.extend_from_slice(q.as_slice());
                }
            }
            last_frame = *run.frames.last().unwrap_or(&last_frame);
            bytes.clear();
            p.pack(&table, &mut bytes);
            ctx.write(&tab_buf, &bytes);
            let gpu_start = Instant::now();
            let mut k0 = 0;
            while k0 < len {
                let k1 = (k0 + per_dispatch).min(len);
                let mut params = Params::new(n_tests as u32);
                (params.a, params.b, params.c) = (k0 as u32, k1 as u32, np as u32);
                ctx.write(&params_buf, bytemuck::bytes_of(&params));
                let mut encoder = ctx.device.create_command_encoder(&Default::default());
                {
                    let mut pass = encoder.begin_compute_pass(&Default::default());
                    pass.set_pipeline(&pipeline);
                    pass.set_bind_group(0, &group, &[]);
                    pass.dispatch_workgroups(groups, 1, 1);
                }
                ctx.queue.submit([encoder.finish()]);
                k0 = k1;
            }
            done += len;

            // The next block on the CPU while the GPU works.
            let next = block.min(n - done);
            if next > 0 {
                let tic = Instant::now();
                record(&mut run, next, t);
                cpu_seconds += tic.elapsed().as_secs_f64();
            }
            let tic = Instant::now();
            ctx.wait();
            wait_seconds += tic.elapsed().as_secs_f64();
            let gpu = gpu_start.elapsed().as_secs_f64();
            if gpu > 0.0 {
                let per_step = gpu / len as f64;
                per_dispatch = ((0.05 / per_step) as usize).clamp(1, block);
            }
            if next == 0 {
                break;
            }
            len = next;
        }

        let st = p.unpack(&ctx.read(&st_buf));
        let flags: Vec<u32> = bytemuck::cast_slice(&ctx.read(&flags_buf)).to_vec();
        run.write_back(states);
        let mut failed = Vec::new();
        for (k, &i) in test_idx.iter().enumerate() {
            if flags[k] != 0 {
                failed.push(i);
                continue;
            }
            let s = &st[6 * k..6 * k + 6];
            let x = Vector3::new(s[0], s[1], s[2]) + last_frame.xc_out;
            let v = Vector3::new(s[3], s[4], s[5]) + last_frame.v_cm_out;
            states[i] = from_pv(&x, &v);
        }
        if !failed.is_empty() {
            // Redo those particles on the CPU with the massive bodies, which take the same steps.
            let keep: Vec<usize> = (0..states.len()).filter(|&i| masses[i] != 0.0).chain(failed.iter().copied()).collect();
            let mut sub: Vec<State> = keep.iter().map(|&i| initial[i]).collect();
            let sub_m: Vec<f64> = keep.iter().map(|&i| masses[i]).collect();
            let central_sub = keep.iter().position(|&i| i == central).unwrap();
            let mut t_sub = 0.0;
            self.cpu.timestep = h;
            self.cpu.steps_blocked(&mut sub, &sub_m, &mut t_sub, central_sub, forces, n);
            for (s, &i) in sub.iter().zip(&keep) {
                if masses[i] == 0.0 {
                    states[i] = *s;
                }
            }
        }

        let mut stats = self.stats.lock().unwrap();
        stats.particle_steps += (n_tests * n) as f64;
        stats.cpu_seconds += cpu_seconds;
        stats.gpu_wait_seconds += wait_seconds;
        stats.total_seconds += start.elapsed().as_secs_f64();
        stats.cpu_fallbacks += failed.len();
    }
}

impl Integrator for GpuWisdomHolman {
    fn step(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>]) {
        self.cpu.timestep = self.timestep;
        self.cpu.step(states, masses, t, forces);
    }

    fn steps(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>], n: usize) {
        let central = central_body(masses);
        let newtonian_only = forces.len() == 1 && forces[0].is_newtonian_gravity();
        let n_tests = masses.iter().filter(|&&m| m == 0.0).count();
        match central {
            Some(central) if newtonian_only && n_tests > 0 && n_tests * n >= self.min_work => self.steps_gpu(states, masses, t, central, forces, n),
            _ => {
                self.cpu.timestep = self.timestep;
                self.cpu.steps(states, masses, t, forces, n);
            }
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

/// One GPU Kepler step for each `(r0, v0, gm, dt)` in `cases`, for testing the kernel's solver
/// against [`universal_kepler_step`](crate::transforms::universal_kepler_step): the position,
/// the velocity, and whether the solver converged.
#[doc(hidden)]
pub fn kepler_steps(ctx: &GpuContext, precision: Precision, cases: &[([f64; 3], [f64; 3], f64, f64)]) -> Vec<([f64; 6], bool)> {
    let src = ctx.source(precision, &kepler_constants(precision), &[include_str!("shaders/kepler.wgsl"), include_str!("shaders/kepler_test.wgsl")]);
    let pipeline = ctx.pipeline(&ctx.module("kepler_test", &src), "main");
    let flat: Vec<f64> = cases.iter().flat_map(|(r, v, gm, dt)| [r[0], r[1], r[2], v[0], v[1], v[2], *gm, *dt]).collect();
    let mut bytes = vec![];
    precision.pack(&flat, &mut bytes);
    let input = ctx.storage(&bytes, false);
    let out = ctx.storage_zeroed((7 * cases.len() * precision.bytes()) as u64);
    let params = ctx.uniform(&Params::new(cases.len() as u32));
    let group = ctx.bind_group(&pipeline, &[&params, &input, &out]);
    ctx.dispatch(&pipeline, &group, (cases.len() as u32).div_ceil(64));
    precision.unpack(&ctx.read(&out)).chunks_exact(7).map(|c| (c[..6].try_into().unwrap(), c[6] == 1.0)).collect()
}
