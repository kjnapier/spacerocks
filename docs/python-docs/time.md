<h1 style="border-bottom: 5px solid white;">Time Module</h1>

### Table of Contents
1. [Overview](#overview)
2. [Primary Classes](#primary-classes)
3. [Methods](#methods)
4. [Examples](#examples)
5. [Notes](#notes)

<h2 style="border-bottom: 3px solid white;">Overview</h2>

The Time module provides functionality for handling astronomical time calculations. It supports multiple time scales (UTC, TDB, TT, TAI) and formats (JD, MJD), with conversion utilities and arithmetic operations.

<h2 style="border-bottom: 3px solid white;">Primary Classes</h2>

### Time
```python
class Time:
    """
    Core class for handling astronomical time calculations.
    
    Example:
        t = Time(2451545.0, "UTC", "JD")
        print(t.calendar())  # "01 Jan 2000"
    """
```

<h2 style="border-bottom: 3px solid white;">Methods</h2>

### Constructor Methods
---

**`new()`**

**Arguments:**
- `epoch`: The epoch value (JD or MJD)
- `timescale`: The timescale ("UTC", "TDB", "TT", "TAI")
- `format`: The time format ("JD", "MJD")

**Returns:**
- New Time instance

*Example:*
```python
time = Time(2451545.0, "UTC", "JD")
```

**`now()`**
```python
@classmethod
def now(cls) -> Time
```
**Returns:**
- Time instance representing the current time in UTC

*Example:*
```python
current_time = Time.now()
```

**`from_fuzzy_str()`**
```python
@classmethod
def from_fuzzy_str(cls, s: str) -> Time
```
**Arguments:**
- `s`: String in format "epoch timescale format" or "now"

**Returns:**
- Time instance parsed from string

*Example:*
```python
time = Time.from_fuzzy_str("2451545.0 UTC JD")
```

**`infer_time_format()`**
```python
@classmethod
def infer_time_format(cls, epoch: float, timescale: Optional[str] = None) -> Time
```
**Arguments:**
- `epoch`: The epoch value
- `timescale`: Optional timescale (defaults to "UTC")

**Returns:**
- Time instance with inferred format (jd or mjd)

*Example:*
```python
time = Time.infer_time_format(2451545.0)
```

### Time Scale Methods
---

**`utc()`**
```python
def utc(self) -> Time
```
**Returns:**
- New Time object in UTC scale

*Example:*
```python
utc_time = time.utc()
```

**`to_utc()`**
```python
def to_utc(self) -> Time
```
Converts time to UTC scale in place.

*Example:*
```python
time.to_utc()
```

Similar methods exist for other timescales:
- **`tdb()`** / **`to_tdb()`**
- **`tt()`** / **`to_tt()`**
- **`tai()`** / **`to_tai()`**

### Format Methods

**`jd()`**
```python
def jd(self) -> float
```
**Returns:**
- Epoch as Julian Date

*Example:*
```python
jd = time.jd()
```

**`mjd()`**
```python
def mjd(self) -> float
```
**Returns:**
- Epoch as Modified Julian Date

*Example:*
```python
mjd = time.mjd()
```

**`calendar()`**
```python
def calendar(self) -> str
```
**Returns:**
- Human-readable calendar date

*Example:*
```python
date = time.calendar()  # "01 Jan 2000"
```

**`iso()`**
```python
def iso(self) -> str
```
**Returns:**
- ISO 8601 formatted string

*Example:*
```python
date = time.iso()  # "2000-01-01T12:00:00.000Z"
```

### Arithmetic Operators
---

Time objects support addition and subtraction of days:
```python
# Addition
tomorrow = time + 1.0

# Subtraction
yesterday = time - 1.0
```

<h2 style="border-bottom: 3px solid white;">Examples</h2>

### Basic Usage
```python
from spacerocks.time import Time

# Create time object
time = Time(2451545.0, "UTC", "JD")

# Convert between scales
tdb_time = time.tdb()
print(f"TDB epoch: {tdb_time.epoch}")

# Format conversions
print(f"Calendar date: {time.calendar()}")
print(f"ISO format: {time.iso()}")

# Informative error messages when using invalid timescale of timeformat 
time_tdb = Time(20424.5, 'tdx', 'mjd')
# Output:
   # ValueError: Invalid timescale: 'tdx'. Did you mean 'tdb'?. Needs to be 'utc', 'tdb', 'tt', or 'tai'.
```

### Time Scale Conversions
```python
# Create time and convert between scales
time = Time(2451545.0, "UTC", "JD")
time.to_tdb()         # Convert to TDB
time.to_tt()          # Convert to TT
time.to_utc()         # Back to UTC

# Or create new objects
tdb = time.tdb()      # New TDB time
tt = time.tt()        # New TT time
```

<h2 style="border-bottom: 3px solid white;">Notes</h2>

- All timescales are case-insensitive
- Calendar format provides human-readable dates
- ISO format follows the 8601 standard
- Invalid inputs raise ValueError with suggestions 
- Method chaining supported for conversions

- UTC before 1972. From 1972, TT − UTC = 32.184 s plus the leap seconds. From 1960 to 1972, UTC
  ran at offset rates from TAI with small steps; spacerocks follows the same table as SOFA/ERFA,
  and applies each step at 0h UTC, when it happened (ERFA spreads it over the preceding day).
  Before 1960 there was no UTC: times are taken as UT, and TT − UT = ΔT from Stephenson, Morrison
  & Hohenkerk (2016, updated 2021; Table S15.2020), about −2 s in 1900, 29 s in 1950 and 33 s in
  1960. `spacerocks.time` agrees with ERFA from 1960 on and with that table before, to the 40 µs
  precision of a Julian date (`validation/time/compare_erfa.py`). SPICE (and so layup, sorcha)
  instead takes 42.184 s for every date before 1972, which is off by 8–44 s in the 20th century.
