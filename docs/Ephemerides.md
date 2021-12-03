# Ephemerides

An `Ephemerides` object is returned from the `observe` method of a `SpaceRock` object. 
Consider you have instantiared a `SpaceRock` object called `rocks`. You can then compute the 
ephemerides like this

```Python
obs = rocks.observe(obscode='W84')
```

The attributes (listed in the following table) are computed lazily in the interest of computational efficiency.

| Observed Parameter                    | Attribute   |
|:--------------------------------------|:------------|
| Right Ascension                       | ra          |
| Declination                           | dec         |
| Apparent Magnitude                    | mag         |
| Ecliptic Latitude                     | b           |
| Ecliptic Longitude                    | l           |
| RA Rate                               | ra_rate     |
| Dec Rate                              | dec_rate    |

## The `hpix` Method

The `hpix` method computed the object's `HEALPix` values, provided `nside` and `nest` parameters.

```Python
hpix64 = obs.hpix(nside=64, nest=True)
```