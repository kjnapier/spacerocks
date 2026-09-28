# Time Module

### Table of Contents
1. [Overview](#overview)
2. [Time Structure](#time-structure)
3. [Constructor Methods](#constructor-methods)
4. [Time Conversions](#time-conversions)
5. [Format Methods](#format-methods)
6. [Arithmetic Operations](#arithmetic-operations)
7. [Notes](#notes)

## Overview

The Time module provides a comprehensive interface for handling astronomical time calculations. It supports multiple time scales (UTC, TDB, TT, TAI) and formats (JD, MJD), with conversion utilities and arithmetic operations.

## Time Structure

### Time Struct
```rust
pub struct Time {
    pub epoch: f64,
    pub timescale: TimeScale,
    pub format: TimeFormat,
}
```

The `Time` struct represents a moment in time with an epoch value in a specific timescale and format.

## Constructor Methods

#### new()
```rust
fn new(epoch: f64, timescale: &str, format: &str) -> Result<Self, TimeError>
```

Creates a new Time instance with specified epoch, timescale, and format.

**Arguments:**
- `epoch`: The epoch value (JD or MJD)
- `timescale`: The timescale ("UTC", "TDB", "TT", "TAI")
- `format`: The time format ("JD", "MJD")

**Returns:**
- `Ok(Time)` if parameters are valid
- `Err(TimeError)` with suggestion if parameters are invalid

**Example:**
```rust
// Create a time in UTC using Julian Date
let t1 = Time::new(2451545.0, "UTC", "JD")?;

// Create a time in TDB using Modified Julian Date
let t2 = Time::new(51544.5, "TDB", "MJD")?;

// Invalid timescale gets helpful error
let t3 = Time::new(2451545.0, "tax", "jd");
// Error: Invalid timescale: 'tax'. Did you mean 'tai'?
```

#### now()
```rust
fn now() -> Self
```

Creates a new Time instance representing the current time.

**Returns:**
- A new `Time` instance in UTC and JD format

**Example:**
```rust
let current_time = Time::now();
println!("Current JD: {}", current_time.jd());
```

#### from_fuzzy_str()
```rust
fn from_fuzzy_str(s: &str) -> Result<Self, TimeError>
```

Creates a Time instance from a string specification.

**Arguments:**
- `s`: String in format "epoch timescale format" or "now"

**Returns:**
- `Ok(Time)` if string is valid
- `Err(TimeError)` if parsing fails

**Example:**
```rust
// Create from explicit string
let t1 = Time::from_fuzzy_str("2451545.0 UTC JD")?;

// Create using current time
let t2 = Time::from_fuzzy_str("now")?;

// Create using MJD
let t3 = Time::from_fuzzy_str("51544.5 TDB MJD")?;
```

#### infer_time_format()
```rust
fn infer_time_format(epoch: f64, timescale: Option<&str>) -> Result<Self, TimeError>
```

Creates a Time instance by inferring the format from the epoch value.

**Arguments:**
- `epoch`: The epoch value
- `timescale`: Optional timescale (defaults to "UTC")

**Returns:**
- `Ok(Time)` with inferred format
- `Err(TimeError)` if parameters are invalid

**Example:**
```rust
// Large epoch infers JD format
let t1 = Time::infer_time_format(2451545.0, Some("UTC"))?;
assert_eq!(t1.format, TimeFormat::JD);

// Small epoch infers MJD format
let t2 = Time::infer_time_format(51544.5, None)?;
assert_eq!(t2.format, TimeFormat::MJD);
```

## Time Conversions

### Creating New Time Objects

#### utc()
```rust
fn utc(&self) -> Time
```

Creates a new Time object converted to UTC timescale.

**Example:**
```rust
let tdb_time = Time::new(2456205.5, "tdb", "jd")?;
let utc_time = tdb_time.utc();
assert!(tdb_time.epoch != utc_time.epoch);  // Epochs differ due to timescale conversion
```

Similar methods exist for other timescales: `tdb()`, `tt()`, and `tai()`.

### In-Place Conversions

#### to_utc()
```rust
fn to_utc(&mut self) -> &mut Self
```

Converts the time object to UTC timescale in place.

**Example:**
```rust
let mut time = Time::new(2456205.5, "tdb", "jd")?;
time.to_utc();
assert_eq!(time.timescale, TimeScale::UTC);

// Chain operations
time.to_utc().to_tdb().to_tt();
```

Similar methods exist for other timescales: `to_tdb()`, `to_tt()`, and `to_tai()`.

## Format Methods

#### jd()
```rust
fn jd(&self) -> f64
```

Returns the epoch as a Julian Date.

**Example:**
```rust
let time = Time::new(51544.5, "UTC", "MJD")?;
let jd = time.jd();  // Returns 2451545.0
```

#### mjd()
```rust
fn mjd(&self) -> f64
```

Returns the epoch as a Modified Julian Date.

**Example:**
```rust
let time = Time::new(2451545.0, "UTC", "JD")?;
let mjd = time.mjd();  // Returns 51544.5
```

#### calendar()
```rust
fn calendar(&self) -> String
```

Returns a human-readable calendar date.

**Example:**
```rust
let time = Time::new(2451545.0, "UTC", "JD")?;
let date = time.calendar();  // Returns "01 Jan 2000"
```

#### iso()
```rust
fn iso(&self) -> String
```

Returns the time in ISO 8601 format.

**Example:**
```rust
let time = Time::new(2451545.0, "UTC", "JD")?;
let iso = time.iso();  // Returns "2000-01-01T12:00:00.000Z"
```

## Arithmetic Operations

### Addition and Subtraction

```rust
// Add days to time
let time = Time::new(2451545.0, "UTC", "JD")?;
let tomorrow = time + 1.0;  // Adds one day
assert_eq!(tomorrow.epoch, 2451546.0);

// Subtract days from time
let yesterday = time - 1.0;  // Subtracts one day
assert_eq!(yesterday.epoch, 2451544.0);

// Add days in place
let mut time = Time::new(2451545.0, "UTC", "JD")?;
time += 1.0;
assert_eq!(time.epoch, 2451546.0);

// Find difference between times
let time1 = Time::new(2451545.0, "UTC", "JD")?;
let time2 = Time::new(2451546.0, "UTC", "JD")?;
let diff: f64 = &time2 - &time1;  // Returns 1.0 (days)
```

## Notes

1. **Timescale Handling**
   - All timescale conversions maintain precision
   - Conversions chain appropriately (e.g., TAI → TT → TDB)
   - Invalid timescale strings receive helpful suggestions

2. **Format Conversions**
   - JD/MJD conversions handle offset automatically
   - Calendar format returns human-readable dates
   - ISO format provides standard timestamp strings

3. **Error Handling**
   - Invalid inputs receive helpful error messages
   - Suggestions provided for close matches
   - Time arithmetic checks for compatible timescales

4. **Time Arithmetic**
   - Addition and subtraction preserve timescale and format
   - Time differences require matching timescales
   - Panic on timescale mismatch for differences

## UTC before 1972

`utc_to_tai`, `tai_to_utc` and everything built on them use `tt_minus_utc(jd_utc)`:

- from 1972: 32.184 s plus the leap seconds (`leapseconds::LEAP_SECONDS`);
- 1960–1972: 32.184 s plus TAI − UTC from UTC's rate-offset table, as SOFA/ERFA's `dat`, with
  each step applied at 0h UTC;
- before 1960 (`UTC_START`), when there was no UTC: the epoch is UT, and TT − UT = ΔT
  (`time::delta_t`, the cubic splines of Stephenson, Morrison & Hohenkerk 2016, Table S15.2020;
  their long-term parabola outside 720 BC–AD 2019).

Results from 1972 on are bit-identical to before. Validated against ERFA and the table in
`validation/time/compare_erfa.py`. SPICE, and so layup and sorcha, take 42.184 s for all dates
before 1972.
