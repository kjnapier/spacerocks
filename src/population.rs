//! Many bodies in one reference plane, about one origin.

use rayon::prelude::*;

use crate::state::{self, Elements, State};
use crate::{Observer, Origin, Properties, ReferencePlane, SpaceRock, Time};
use crate::observing::{apparent, Apparent};

type BoxError = Box<dyn std::error::Error + Send + Sync>;

/// Many bodies stored as a dense state array plus metadata.
///
/// The reference plane and origin are shared by every body, so they are stored once and checked
/// when bodies are added. Epochs (TDB Julian dates), names and physical properties are per body
/// and kept in separate columns; names and properties are never touched by the numerical code. The physics is done by the
/// free functions in [`crate::state`], mapped over [`Population::states`].
///
/// [`SpaceRock`] stays the one-body type: [`Population::push`] takes one and
/// [`Population::get`] builds one on demand.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Population {
    pub reference_plane: ReferencePlane,
    pub origin: Origin,
    /// `[x, y, z, vx, vy, vz]` per body, AU and AU/day.
    pub states: Vec<State>,
    /// TDB Julian date of each state.
    pub epochs: Vec<f64>,
    pub names: Vec<String>,
    pub properties: Vec<Option<Properties>>,
}

impl Population {
    /// An empty population in `reference_plane` about `origin`.
    pub fn new(reference_plane: ReferencePlane, origin: Origin) -> Population {
        Population { reference_plane, origin, ..Default::default() }
    }

    /// Collect rocks into a population that takes the first rock's reference plane and origin
    /// (J2000 and SSB if there are none); see [`Population::push`] for how the others are
    /// checked.
    pub fn from_rocks<I: IntoIterator<Item = SpaceRock>>(rocks: I) -> Result<Population, BoxError> {
        let mut rocks = rocks.into_iter().peekable();
        let mut pop = match rocks.peek() {
            Some(r) => Population::new(r.reference_plane.clone(), r.origin.clone()),
            None => Population::default(),
        };
        for rock in rocks {
            pop.push(rock)?;
        }
        Ok(pop)
    }

    pub fn len(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// Add a rock. A rock in another reference plane is rotated into the population's (an exact, local
    /// transformation), and a rock about another origin is an error, since changing origin
    /// needs an ephemeris: move it first with [`SpaceRock::to_ssb`] or [`SpaceRock::to_helio`].
    pub fn push(&mut self, mut rock: SpaceRock) -> Result<(), BoxError> {
        if rock.origin != self.origin {
            return Err(format!(
                "{} has origin {} but the collection's origin is {}; change its origin before adding it",
                rock.name, rock.origin, self.origin
            )
            .into());
        }
        if rock.reference_plane != self.reference_plane {
            rock.change_reference_plane(self.reference_plane.as_str()).map_err(|e| e.to_string())?;
        }
        self.states.push(rock.state());
        self.epochs.push(rock.epoch.tdb().jd());
        self.names.push(rock.name);
        self.properties.push(rock.properties);
        Ok(())
    }

    /// The `i`-th body as a [`SpaceRock`].
    pub fn get(&self, i: usize) -> Option<SpaceRock> {
        let (position, velocity) = state::pv(self.states.get(i)?);
        Some(SpaceRock {
            name: self.names[i].clone(),
            epoch: tdb(self.epochs[i]),
            reference_plane: self.reference_plane.clone(),
            origin: self.origin.clone(),
            position,
            velocity,
            properties: self.properties[i].clone(),
        })
    }

    /// Index of the first body called `name`.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }

    /// Every body as a [`SpaceRock`].
    pub fn to_rocks(&self) -> Vec<SpaceRock> {
        (0..self.len()).map(|i| self.get(i).unwrap()).collect()
    }

    /// Replace the states and epochs with those of `rocks` (same order and length, e.g. after
    /// moving them with [`crate::batch`]). The rocks must still be in this population's
    /// reference plane and about its origin.
    pub fn update_from_rocks(&mut self, rocks: &[SpaceRock]) -> Result<(), BoxError> {
        if rocks.len() != self.len() {
            return Err("number of rocks does not match the population".into());
        }
        for (i, rock) in rocks.iter().enumerate() {
            if rock.reference_plane != self.reference_plane || rock.origin != self.origin {
                return Err(format!("{} changed reference plane or origin", rock.name).into());
            }
            self.states[i] = rock.state();
            self.epochs[i] = rock.epoch.tdb().jd();
        }
        Ok(())
    }

    /// The bodies where `mask` is true.
    pub fn filter(&self, mask: &[bool]) -> Result<Population, BoxError> {
        if mask.len() != self.len() {
            return Err("Mask length must match the number of rocks.".into());
        }
        let idx: Vec<usize> = (0..self.len()).filter(|&i| mask[i]).collect();
        Ok(self.select(&idx))
    }

    /// The bodies at `idx`, in that order, with the same plane and origin.
    pub fn select(&self, idx: &[usize]) -> Population {
        Population {
            reference_plane: self.reference_plane.clone(),
            origin: self.origin.clone(),
            states: idx.iter().map(|&i| self.states[i]).collect(),
            epochs: idx.iter().map(|&i| self.epochs[i]).collect(),
            names: idx.iter().map(|&i| self.names[i].clone()).collect(),
            properties: idx.iter().map(|&i| self.properties[i].clone()).collect(),
        }
    }

    /// Osculating elements of every body about the origin, in one pass per body.
    pub fn elements(&self) -> Vec<Elements> {
        let mu = self.origin.mu();
        self.states.par_iter().map(|s| state::elements(s, mu)).collect()
    }

    /// Apply `f` to every state with the origin's μ, in parallel. The building block for
    /// column getters like the semi-major axis.
    pub fn map<T: Send, F: Fn(&State, f64) -> T + Sync + Send>(&self, f: F) -> Vec<T> {
        let mu = self.origin.mu();
        self.states.par_iter().map(|s| f(s, mu)).collect()
    }

    /// Rotate every body into `reference_plane`.
    pub fn change_reference_plane(&mut self, reference_plane: &ReferencePlane) -> Result<(), BoxError> {
        if *reference_plane == self.reference_plane {
            return Ok(());
        }
        let rot = state::rotation(&self.reference_plane, reference_plane).map_err(|e| e.to_string())?;
        self.states.par_iter_mut().for_each(|s| *s = state::rotate(s, &rot));
        self.reference_plane = reference_plane.clone();
        Ok(())
    }

    /// Move every body along its two-body orbit about the origin to `epoch`.
    pub fn analytic_propagate(&mut self, epoch: &Time) -> Result<(), BoxError> {
        let mu = self.origin.mu();
        let t = epoch.tdb().jd();
        let names = &self.names;
        self.states
            .par_iter_mut()
            .zip(self.epochs.par_iter())
            .enumerate()
            .try_for_each(|(i, (s, ep))| {
                *s = state::kepler_step(s, mu, t - ep).map_err(|e| format!("{}: {}", names[i], e))?;
                Ok::<(), String>(())
            })?;
        self.epochs.iter_mut().for_each(|ep| *ep = t);
        Ok(())
    }

    /// Observable quantities of every body as seen by `observer` (see [`SpaceRock::apparent`]
    /// for the conventions). The observer must be in the population's reference plane, have a
    /// velocity, and be at every body's epoch to within 1 µs.
    pub fn apparent(&self, observer: &Observer) -> Result<Vec<Apparent>, BoxError> {
        let obs_vel = observer.velocity.ok_or("Observer velocity is required to compute rates")?;
        if observer.reference_plane != self.reference_plane {
            return Err("Observer and rocks have different reference planes".into());
        }
        let t_obs = observer.epoch.tdb().jd();
        let sun = observer.sun();
        (0..self.len())
            .into_par_iter()
            .map(|i| {
                if (self.epochs[i] - t_obs).abs() > 1e-6 / 86400.0 {
                    return Err(format!("{}: Observer and SpaceRock have different epochs", self.names[i]).into());
                }
                let (h, g) = match &self.properties[i] {
                    Some(p) => (p.absolute_magnitude, p.gslope.unwrap_or(0.15)),
                    None => (None, 0.15),
                };
                let (r, v) = state::pv(&self.states[i]);
                Ok(apparent(&r, &v, &observer.position, &obs_vel, &sun, h, g))
            })
            .collect()
    }
}

/// A TDB Julian date as a [`Time`].
pub(crate) fn tdb(jd: f64) -> Time {
    Time::new(jd, "tdb", "jd").expect("TDB JD is a valid time")
}
