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
| right ascension                       | ra          |
| declination                           | dec         |
| apparent magnitude                    | mag         |
| ecliptic latitude                     | b           |
| ecliptic longitude                    | l           |
| ra rate                               | ra_rate     |
| dec rate                              | dec_rate    |

| Physical Parameter                               | Attribute   |
|:-------------------------------------------------|:------------|
| equatorial x position with respect to observer   | x           |
| equatorial y position with respect to observer   | y           |
| equatorial z position with respect to observer   | z           |
| equatorial vx velocity with respect to observer  | vx          |
| equatorial vy velocity with respect to observer  | vy          |
| equatorial vz velocity with respect to observer  | vz          |
| distance from observer                           | delta       |

| Physical Property                     | Attribute   |
|:--------------------------------------|:------------|
| absolute magnitude                    | H           |
| phase slope constant                  | G           |


| Metadata                              | Attribute   |
|:--------------------------------------|:------------|
| name                                  | name        |
| epoch                                 | epoch       |

## The `hpix` Method

The `hpix` method computes the object's `HEALPix` values, provided `nside` and `nest` parameters.

```Python
hpix64 = obs.hpix(nside=64, nest=True)
```