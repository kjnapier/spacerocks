import pandas as pd
from astropy.coordinates import earth
import numpy as np
pd.options.mode.chained_assignment = None

observatories = pd.read_json('/Users/kjnapier/research/spacerocks/spacerocks/obscodes_extended.json').T
observatories['obscode'] = observatories.index
todrop = observatories[observatories.Latitude.isnull().values].index
observatories = observatories.drop(todrop)
observatories['Latitude'] = np.degrees(observatories.Latitude.values.astype(float))
observatories['Longitude'][observatories.Longitude.values > 180] -= 360
observatories = observatories.rename(columns={'Latitude': 'lat',
                                              'Longitude': 'lon',
                                              'Geocentric_dist': 'geocentric_dist',
                                              'Name': 'name'})
observatories = observatories.drop(columns=['cos', 'sin'])

elevations = np.zeros(len(observatories))
# This function isn't vectorized???
for idx in range(len(observatories)):
    x, y, z = earth.EarthLocation(lat=observatories.lat.values[idx],
                                  lon=observatories.lon.values[idx],
                                  height=0).value
    ellipsoid_edge = np.sqrt(x**2 + y**2 + z**2)
    elevations[idx] = 6378137 * observatories.geocentric_dist.values[idx] - ellipsoid_edge

observatories['elevation'] = elevations
observatories = observatories.drop(columns=['geocentric_dist'])
observatories = observatories[['obscode', 'lat', 'lon', 'elevation', 'name']]
observatories['lat'] = [str(abs(ii)) + ' N' if ii > 0 else str(abs(ii)) + ' S' \
                        for ii in observatories.lat]
observatories['lon'] = [str(abs(ii)) + ' E' if ii > 0 else str(abs(ii)) + ' W' \
                        for ii in observatories.lon]
observatories.to_csv('/Users/kjnapier/research/spacerocks/spacerocks/observatories.csv', index=False)
