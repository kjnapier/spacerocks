# Orbital Transform Functions

### Table of Contents
1. [Overview](#overview)
2. [Orbit Types](#orbit-types)
3. [Anomaly Transforms](#anomaly-transforms)
4. [State Vector Operations](#state-vector-operations)
5. [Examples](#examples)
6. [Notes](#notes)

## Overview

The `transforms` module provides functions for orbital mechanics calculations, focusing on conversions between different types of orbital anomalies and state vectors. All angular units are in radians, distances in astronomical units (AU), and time in days.

## Orbit Types

```rust
pub enum OrbitType {
    Hyperbolic,   // e > 1
    Parabolic,    // e ≈ 1
    Elliptical,   // 0 < e < 1
    Circular,     // e ≈ 0
    Radial        // Special case
}
```

The module also defines an error type for handling invalid inputs:

```rust
pub enum OrbitError {
    NegativeEccentricity(f64),
}
```

#### from_eccentricity()
```rust
fn from_eccentricity(e: f64, threshold: f64) -> Result<OrbitType, OrbitError>
```

Determines orbit type based on eccentricity value. Uses a threshold parameter to classify near-circular and near-parabolic orbits.

**Arguments:**
- `e`: Orbital eccentricity
- `threshold`: Tolerance for determining near-circular or near-parabolic orbits (typically 1e-10)

**Returns:**
- `Ok(OrbitType)` with the appropriate orbit classification
- `Err(OrbitError::NegativeEccentricity)` if e < 0

**Example:**
```rust
let orbit_type = OrbitType::from_eccentricity(0.5, 1e-10)?;
```

## Anomaly Transforms

#### calc_conic_anomaly_from_mean_anomaly()
```rust
fn calc_conic_anomaly_from_mean_anomaly(
    e: f64,
    mean_anomaly: f64
) -> Result<f64, OrbitError>
```

Converts mean anomaly to conic anomaly using orbit-type-specific methods:
- Circular orbits: Returns mean anomaly directly
- Elliptical orbits: Uses Newton-Raphson iteration with third-order corrections
- Parabolic orbits: Uses analytical solution
- Hyperbolic orbits: Uses iterative solution with convergence check

The elliptical case includes special handling for mean anomalies greater than π.

**Arguments:**
- `e`: Orbital eccentricity
- `mean_anomaly`: Mean anomaly in radians

**Example:**
```rust
let eccentric_anomaly = calc_conic_anomaly_from_mean_anomaly(0.5, 0.5)?;
```

#### calc_conic_anomaly_from_true_anomaly()
```rust
fn calc_conic_anomaly_from_true_anomaly(
    e: f64,
    true_anomaly: f64
) -> Result<f64, Box<dyn std::error::Error>>
```

Converts true anomaly to conic anomaly using specific formulas for each orbit type:
- Circular orbits: Returns true anomaly directly
- Elliptical orbits: Uses arctan formula with eccentricity corrections
- Parabolic orbits: Uses tangent half-angle formula
- Hyperbolic orbits: Uses hyperbolic arctangent formula

**Arguments:**
- `e`: Orbital eccentricity
- `true_anomaly`: True anomaly in radians

**Example:**
```rust
let eccentric_anomaly = calc_conic_anomaly_from_true_anomaly(0.3, 1.2)?;
```

#### calc_mean_anomaly_from_conic_anomaly()
```rust
fn calc_mean_anomaly_from_conic_anomaly(
    e: f64,
    conic_anomaly: f64
) -> Result<f64, OrbitError>
```

Converts conic anomaly to mean anomaly using Kepler's equations:
- Circular orbits: M = E
- Elliptical orbits: M = E - e sin(E)
- Parabolic orbits: M = B - B³/3
- Hyperbolic orbits: M = e sinh(H) - H

Where E is eccentric anomaly, B is parabolic eccentric anomaly, and H is hyperbolic eccentric anomaly.

**Arguments:**
- `e`: Orbital eccentricity
- `conic_anomaly`: Conic anomaly in radians

**Example:**
```rust
let mean_anomaly = calc_mean_anomaly_from_conic_anomaly(0.5, 0.5)?;
```

## State Vector Operations

Orbital elements from a state are computed by `spacerocks::state::elements(&state, mu)`
(all elements in one pass) and the single-element functions next to it
(`state::eccentricity`, `state::inclination`, ...).

## Examples

### Complete Anomaly Conversion
```rust
use spacerocks::transforms;

// Convert from mean anomaly to true anomaly for an elliptical orbit
let e = 0.7;
let mean_anomaly = 0.8;

// First get eccentric anomaly
let eccentric_anomaly = transforms::calc_conic_anomaly_from_mean_anomaly(e, mean_anomaly)?;

// Then convert to true anomaly
let true_anomaly = transforms::calc_conic_anomaly_from_true_anomaly(e, eccentric_anomaly)?;
```

### Error Handling
```rust
use spacerocks::transforms;

// Handle invalid eccentricity
let e = -0.1;
match transforms::calc_conic_anomaly_from_mean_anomaly(e, 0.0) {
    Ok(_) => println!("Valid calculation"),
    Err(e) => println!("Error: {}", e),
}
```

## Notes

1. **Orbit Classification**
   - Eccentricity determines orbit type (circular, elliptical, parabolic, hyperbolic)
   - Small threshold value (typically 1e-10) used for circular/parabolic classification
   - Negative eccentricities are invalid and return errors

2. **Numerical Methods**
   - Newton-Raphson iteration used for elliptical orbits with third-order corrections
   - Convergence tolerance of 1e-15 for iterative solutions
   - Special analytical solutions for circular and parabolic cases
   - Hyperbolic orbits use specialized iteration methods

3. **Units and Conventions**
   - All angles in radians
   - Distances in astronomical units (AU)
   - Times in days
   - Velocities in AU/day
   - State vectors use nalgebra Vector3 type