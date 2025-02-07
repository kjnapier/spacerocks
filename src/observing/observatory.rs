use crate::SpaceRock;

use crate::observing::observer::Observer;
use crate::constants::{DEG_TO_RAD, M_TO_AU, EQUAT_RAD};
use crate::OBSERVATORIES;
use crate::time::Time;

use nalgebra::Vector3;


/// Represents different types of astronomical observatories
#[derive(Clone, Debug, PartialEq)]
pub enum Observatory {
     /// Fixed ground-based observatory with known position on Earth
    GroundObservatory { obscode: String, lon: f64, lat: f64, rho: f64 },
     /// Space-based telescope with position from SPICE kernels
    SpaceTelecope { name: String },
    /// Observing from a SpaceRock
    SpaceRockObservatory { rock: SpaceRock }
}

impl Observatory {

    /// Create a new Observatory from an observatory code.
    ///
    /// # Arguments
    ///
    /// * `obscode` - The observatory code.
    ///
    /// # Returns
    ///
    /// * `Observatory` - The Observatory object.
    pub fn from_obscode(obscode: &str) -> Result<Self, &'static str> {
        let obscode = obscode.to_uppercase();
        match OBSERVATORIES.get(&obscode) {
            Some(obs) => {
                let lon = obs.0;
                let lat = obs.2.atan2(obs.1);
                let rho = (obs.1 * obs.1 + obs.2 * obs.2).sqrt();
                return Ok(Observatory::GroundObservatory { obscode: obscode, lon, lat, rho })
            },
            None => {
                return Err("Observatory not found")
            }
        }
    }

    /// Create a new Observatory from a name. 
    /// The name should usually the name of a space telescope, but it can be anything that 
    /// is loaded into the SPICE kernel.
    ///
    /// # Arguments
    ///
    /// * `name` - Name of the observatory.
    ///
    /// # Returns
    ///
    /// * `Observatory` - The Observatory object.
    pub fn from_name(name: &str) -> Self {
        Observatory::SpaceTelecope { name: name.to_string() }
    }
   
    /// Get the Observer at a specific time.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The time to get the observer at.
    ///
    /// # Returns
    ///
    /// * `Result<Observer, Box<dyn std::error::Error>>` - The Observer object.
    pub fn at(&self, epoch: &Time, reference_plane: &str, origin: &str) -> Result<Observer, Box<dyn std::error::Error>> {
        match self {
            Observatory::GroundObservatory { obscode: _, lon, lat, rho } => {
                let mut earth = SpaceRock::from_spice("earth", epoch, reference_plane, origin)?;
                let rho_sin_lat = lat.sin() * rho;
                let rho_cos_lat = lat.cos() * rho;
                
                let delta_et = 10.0;
                let et = spice::str2et(&format!("JD{epoch} UTC", epoch=epoch.utc().jd()));     
                let m: nalgebra::Matrix3<f64> = spice::pxform("ITRF93", reference_plane, et).into();
                let mp: nalgebra::Matrix3<f64> = spice::pxform("ITRF93", reference_plane, et + delta_et).into();
                let mm: nalgebra::Matrix3<f64> = spice::pxform("ITRF93", reference_plane, et - delta_et).into();
                // transpose the matrix
                let m = m.transpose();
                let mp = mp.transpose();
                let mm = mm.transpose();

                let ox = rho_cos_lat * lon.cos();
                let oy = rho_cos_lat * lon.sin();
                let oz = rho_sin_lat;
                let obsVec = Vector3::new(ox, oy, oz);

                // println!("{:?}", obsVec);

                let mVec = m * obsVec * EQUAT_RAD * M_TO_AU;
                let mVecp = mp * obsVec * EQUAT_RAD * M_TO_AU;
                let mVecm = mm * obsVec * EQUAT_RAD * M_TO_AU;
                let d_vel = (mVecp - mVecm) / (2.0 * delta_et / 86400.0);
                earth.position += mVec;
                earth.velocity += d_vel;
                
                // let [d_pos, d_vel] = compute_topocentric_correction(*lon, rho_sin_lat, rho_cos_lat, epoch.jd());
                // earth.position += d_pos;
                // earth.velocity += d_vel;
                Ok(Observer { spacerock: earth, observatory: self.clone() })
            },
            Observatory::SpaceTelecope { name } => {
                let rock = SpaceRock::from_spice(&name, epoch, reference_plane, origin)?;
                Ok(Observer { spacerock: rock, observatory: self.clone() })
            }
            _ => {
                return Err("at not implemented for this observatory type".into())
            }
        }
    }

    /// Get the name of the Observatory.
    ///
    /// # Returns
    ///
    /// * `String` - The name of the Observatory.
    pub fn name(&self) -> String {
        match self {
            Observatory::GroundObservatory { obscode, .. } => obscode.clone(),
            Observatory::SpaceTelecope { name } => name.clone(),
            Observatory::SpaceRockObservatory { rock } => rock.name.clone()
        }
    }

    /// Get the longitude of the Observatory.
    ///
    /// # Returns
    ///
    /// * `Option<f64>` - The longitude of the Observatory.
    pub fn lon(&self) -> Option<f64> {
        match self {
            Observatory::GroundObservatory { lon, .. } => Some(*lon),
            Observatory::SpaceTelecope { .. } => None,
            Observatory::SpaceRockObservatory { .. } => None,
        }
    }

    /// Get the latitude of the Observatory.
    ///
    /// # Returns
    ///
    /// * `Option<f64>` - The latitude of the Observatory.
    pub fn lat(&self) -> Option<f64> {
        match self {
            Observatory::GroundObservatory { lat, .. } => Some(*lat),
            Observatory::SpaceTelecope { .. } => None,
            Observatory::SpaceRockObservatory { .. } => None,
        }
    }

    /// Get the distance of the observatory from the Geocenter.
    ///
    /// # Returns
    ///
    /// * `Option<f64>` - The rho of the Observatory.
    pub fn rho(&self) -> Option<f64> {
        match self {
            Observatory::GroundObservatory { rho, .. } => Some(*rho),
            Observatory::SpaceTelecope { .. } => None,
            Observatory::SpaceRockObservatory { .. } => None,
        }
    }

}

/// Computes local sidereal time for a given epoch and longitude
/// 
/// # Arguments
/// * `epoch` - Julian date
/// * `lon` - Longitude in radians
fn compute_local_sidereal_time(epoch: f64, lon: f64) -> f64 {
    let t = (epoch - 2451545.0) / 36525.0;
    let mut theta = 280.46061837 + 360.98564736629 * (epoch - 2451545.0) + (0.000387933 * t * t) - (t * t * t / 38710000.0);
    theta *= DEG_TO_RAD;
    return theta + lon
}

/// Computes sidereal rotation rate at given epoch
/// 
/// # Arguments
/// * `epoch` - Julian date
fn sidereal_rate(epoch: f64) -> f64 {
    let t = (epoch - 2451545.0) / 36525.0;
    let tprime = 1.0 / 36525.0;
    let theta_dot = 360.98564736629 + 2.0 * 0.000387933 * t * tprime + 3.0 * t * t * tprime / 38710000.0;
    return theta_dot * DEG_TO_RAD
}

/// Calculates position and velocity corrections for Earth-based observers
/// 
/// # Arguments
/// * `lon` - Observer longitude in radians
/// * `rho_sin_lat` - ρsin(φ) where ρ is distance from Earth center
/// * `rho_cos_lat` - ρcos(φ) where ρ is distance from Earth center
/// * `epoch` - Julian date
/// 
/// # Returns
/// * Position correction [dx, dy, dz] in AU
/// * Velocity correction [vx, vy, vz] in AU/day
fn compute_topocentric_correction(lon: f64, rho_sin_lat: f64, rho_cos_lat: f64, epoch: f64) -> [Vector3<f64>; 2] {

    let phi = compute_local_sidereal_time(epoch, lon);
    let phi_rate = sidereal_rate(epoch);
    
    let sin_lon = phi.sin();
    let cos_lon = phi.cos();

    let sin_lon_prime = cos_lon * phi_rate;
    let cos_lon_prime = -sin_lon * phi_rate;
        
    let dx = rho_cos_lat * cos_lon * EQUAT_RAD;
    let dy = rho_cos_lat * sin_lon * EQUAT_RAD;
    let dz = rho_sin_lat * EQUAT_RAD;

    let dvx = cos_lon_prime * rho_cos_lat * EQUAT_RAD;
    let dvy = sin_lon_prime * rho_cos_lat * EQUAT_RAD;
    let dvz = 0.0;

    let d_pos = Vector3::new(dx, dy, dz) * M_TO_AU; // AU
    let d_vel = Vector3::new(dvx, dvy, dvz) * M_TO_AU; // AU/day

    return [d_pos, d_vel];

}


// earth_latest_high_prec.bpc
// pxform

// def barycentricObservatoryRates(et, obsCode, observatories, Rearth=RADIUS_EARTH_KM, delta_et=10):
//     """
//     Computes the position and rate of motion for the observatory in barycentric coordinates

//     Parameters
//     ----------
//     et: float
//         JPL ephemeris time
//     obsCode: str
//         MPC observatory code
//     observatories: Observatory
//         Observatory object with spherical representations for the obsCode
//     Rearth: float
//         Radius of the Earth (default is RADIUS_EARTH_KM)
//     delta_et: float
//         Difference in ephemeris time (in days) to derive the rotation matrix from the fixed Earth equatorial frame to J2000 (default: 10)
//     Returns
//     -------
//      : array
//         Position of the observatory (baricentric)
//      : array
//         Velocity of the observatory (baricentric)
//     """
//     # This JPL's quoted Earth radius (km)
//     # et is JPL's internal time
//     # Get the barycentric position of Earth
//     posvel, _ = spice.spkezr("EARTH", et, "J2000", "NONE", "SSB")
//     pos = posvel[0:3]
//     vel = posvel[3:6]
//     # Get the matrix that rotates from the Earth's equatorial body fixed frame to the J2000 equatorial frame.
//     m = spice.pxform("ITRF93", "J2000", et)
//     mp = spice.pxform("ITRF93", "J2000", et + delta_et)
//     mm = spice.pxform("ITRF93", "J2000", et - delta_et)
//     # Get the MPC's unit vector from the geocenter to
//     # the observatory
//     obsVec = observatories.ObservatoryXYZ[obsCode]
//     obsVec = np.array(obsVec)
//     # Carry out the rotation and scale
//     mVec = np.dot(m, obsVec) * Rearth
//     mVecp = np.dot(mp, obsVec) * Rearth
//     mVecm = np.dot(mm, obsVec) * Rearth
//     return pos + mVec, vel + (mVecp - mVecm) / (2 * delta_et)