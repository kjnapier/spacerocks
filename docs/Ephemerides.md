# Ephemerides

An `Ephemerides` object is returned from the `observe` method of a `SpaceRock` object. 
Consider you have instantiared a `SpaceRock` object called `rocks`. You can then compute the 
ephemerides like this

```Python
obs = rocks.observe(obscode='W84')
```

The attributes (listed in the following tables) are computed lazily in the interest of computational efficiency.

| Observed Parameter                    | Attribute   |
|:--------------------------------------|:------------|
| Right Ascension                       | ra          |
| Declination                           | dec         |
| Apparent Magnitude                    | mag         |
| Ecliptic Latitude                     | b           |
| Ecliptic Longitude                    | l           |
| RA Rate                               | ra_rate     |
| Dec Rate                              | dec_rate    |

| Physical Parameter                    | Attribute   |
|:--------------------------------------|:------------|
| x position with respect to observer   | x           |
| y position with respect to observer   | y           |
| z position with respect to observer   | z           |
| vx velocity with respect to observer  | vx          |
| vy velocity with respect to observer  | vy          |
| vz velocity with respect to observer  | vz          |

| Physical Property                     | Attribute   |
|:--------------------------------------|:------------|
| Absolute Magnitude                    | H           |

| Metadata                              | Attribute   |
|:--------------------------------------|:------------|
| Name                                  | name        |

## The `hpix` Method

The `hpix` method computes the object's `HEALPix` values, provided `nside` and `nest` parameters.

```Python
hpix64 = obs.hpix(nside=64, nest=True)
```