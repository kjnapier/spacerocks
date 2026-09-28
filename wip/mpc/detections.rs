use reqwest::Error;
use serde::Deserialize;
use std::env;

use crate::observing::{Observation, Observatory};
use crate::time::Time;

const MPC_API: &str = "https://data.minorplanetcenter.net/api/get-obs";

// Only the fields we care about.
#[derive(Debug, Deserialize)]
struct MPCObservation {
    #[serde(deserialize_with = "string_to_f64_option")]
    ra: Option<f64>,
    #[serde(deserialize_with = "string_to_f64_option")]
    dec: Option<f64>,
    #[serde(deserialize_with = "string_to_f64_option")]
    mag: Option<f64>,
    stn: Option<String>,
    obstime: Option<String>,
    #[serde(deserialize_with = "string_to_f64_option")]
    pos1: Option<f64>,
    #[serde(deserialize_with = "string_to_f64_option")]
    pos2: Option<f64>,
    #[serde(deserialize_with = "string_to_f64_option")]
    pos3: Option<f64>,
}


// The top-level response from the API: the first element contains our data.
#[derive(Debug, Deserialize)]
struct GetObsResponse {
    ADES_DF: Option<Vec<MPCObservation>>,
}

async pub fn fetch_detections(designation: &str) -> Result<reqwest::Response, Error> {
    let binding = "ADES_DF".to_string();

    let mut params = std::collections::HashMap::new();
    params.insert("desigs", vec![designation]);
    params.insert("output_format", vec![&binding]);

    let client = reqwest::Client::new();
    let response = client.get(MPC_API).json(&params).send().await?;

    Ok(response)
}

pub async fn get_detections(designation: &str) -> Result<Vec<Observation>, Error> {
    let response = fetch_detections(designation).await?;
    let mut observations = Vec<Observation>::new();
    
    if response.status().is_success() {
        let api_response: (GetObsResponse, i64) = response.json().await?;
        if let Some(observations) = api_response.0.ADES_DF {
            for (i, obs) in observations.iter().enumerate() {
                let ra = obs.ra.unwrap();
                let dec = obs.dec.unwrap();
                let mag = obs.mag.unwrap_or(None);
                let stn = obs.stn.as_ref().unwrap().to_string();
                let obstime = obs.obstime.as_ref().unwrap().to_string();

                let epoch = Time::from_isot(obstime);
                let observatory = Observatory::from_obscode(&stn);
                let observer = observatory.at(epoch);
                let observation = Observation::from_astrometry(epoch, ra, dec, observer, None, mag, None);
                observations.push(observation);
            }
        } else {
            println!("No observations found.");
        }
    } else {
        eprintln!("Request failed with status: {}", response.status());
    }

    Ok(observations)
}



fn string_to_f64_option<'de, D>(deserializer: D) -> Result<Option<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    let opt = Option::<String>::deserialize(deserializer)?;
    if let Some(s) = opt {
        s.parse::<f64>().map(Some).map_err(serde::de::Error::custom)
    } else {
        Ok(None)
    }
}