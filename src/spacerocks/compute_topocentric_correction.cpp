#include "spacerocks.hpp"

struct Vector3 compute_topocentric_correction(double earth_lat, double earth_lon, double elevation, double epoch) {

    Vector3 vec;

    double sin_lat, sin_lon, cos_lat, cos_lon, denom, C_geo, S_geo, lon;
        
    sin_lat = sin(earth_lat);
    cos_lat = cos(earth_lat);
    lon = compute_lst(epoch, earth_lon);

    sin_lon = sin(lon);
    cos_lon = cos(lon);

    denom = O_M_FLATTEN * sin_lat;
    denom = cos_lat * cos_lat + denom*denom;
    C_geo = 1 / sqrt(denom);
    S_geo = O_M_FLATTEN * O_M_FLATTEN * C_geo;
    C_geo = C_geo * EQUAT_RAD + elevation;
    S_geo = S_geo * EQUAT_RAD + elevation;

    vec.x = C_geo * cos_lat * cos_lon;
    vec.y = C_geo * cos_lat * sin_lon;
    vec.z = S_geo * sin_lat;

    return vec;

}