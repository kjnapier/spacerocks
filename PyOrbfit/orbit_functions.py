from skyfield.api import Topos, Loader
load = Loader('./Skyfield-Data', expire=False, verbose=False)
planets = load('de423.bsp')
earth = planets['earth']
from astropy.constants import c
c = c.to(u.au / u.year)
GDEGENERATE = 0.2

ts = load.timescale()

# Returns xyz position of ground-based observatory.
# Need to implement obscode! Hard-coded for DECam
def earth3d(self):
    t = ts.tt(jd=self.obsdate)
    earth += Topos('30.169 S', '70.804 W', elevation_m=2200)
    return earth.at(t).position.au * u.au

def kbo2d_linear(self):
    x_earth, y_earth, z_earth = earth3d()

    t_emit = self.obsdate - z_earth / c

    X = self.a + self.adot * t_emit - self.g * x_earth - self.gdot * (self.adot * t_emit**2 - self.g * x_earth * t_emit)
    Y = self.b + self.bdot * t_emit - self.g * y_earth - self.gdot * (self.bdot * t_emit**2 - self.g * y_earth * t_emit)

    dX = np.array([np.ones(nobs), t_emit, np.zeros(nobs), np.zeros(nobs),
                   -x_earth, self.g * x_earth * t_emit - self.adot * t_emit**2])
    dY = np.array([np.zeros(nobs), np.zeros(nobs), np.ones(nobs), t_emit,
                   -y_earth, self.g * y_earth * t_emit - self.bdot * t_emit**2])

    return X, Y, dX, dY

def equatorial_to_tangent(self):
    β = np.arcsin(np.cos(ε) * np.sin(self.dec) -  \
        np.sin(ε) * np.cos(self.dec) * np.sin(self.ra))
    λ = np.arctan((np.cos(ε) * np.cos(self.dec) * np.sin(self.ra) + \
        np.sin(ε) * np.sin(self.dec)) / (np.cos(self.dec) * np.cos(self.ra)))

    θx = (np.cos(β) * np.sin(λ - λ[0])) / (np.sin(b[0]) * np.sin(β) - \
              np.cos(β[0]) * np.cos(β) * np.cos(λ - λ[0]))
    θy = (np.cos(β[0]) * np.cos(β) - np.sin(β[0]) * np.sin(β) * np.cos(λ - λ[0]))/ \
              (np.sin(β[0]) * np.sin(β) + np.cos(β[0]) * np.cos(β) * np.cos(λ - λ[0]))

    dβ_ddec = (np.cos(self.dec)*np.cos(ε) + np.sin(self.ra)*np.sin(self.dec)*np.sin(ε))**2 / \
              (1 - (np.cos(ε)*np.sin(self.dec) - np.cos(self.dec)*np.sin(self.ra)*np.sin(ε))**2)

    dβ_dra = (np.cos(self.ra) * np.cos(self.dec) * np.sin(ε))**2 / \
             (1 - (np.cos(ε)*np.sin(self.dec) - np.cos(self.dec)*np.sin(self.ra)*np.sin(ε))**2)

    return θx, θy, σθx, σθy

def prelim_fit(self):

    X, Y, dX, dY = kbo2d_linear()

    beta = np.dot(self.θx * self.σθx**2, dX) + np.dot(self.θy * self.σθy**2, dY)
    alpha = np.matmul(dX.T, dX * np.tile(1/self.σθx**2, (6, 1)).T) + \
            np.matmul(dY.T, dY * np.tile(1/self.σθy**2, (6, 1)).T)

    try:
        cov = np.linalg.inv(alpha)
    except:
        raise ValueError('The covariance matrix is singular.')

    a, adot, b, bdot, g, _ = np.matmul(alpha, beta)
    gdot = 0

    cov[5], cov.T[5] = np.zeros(6), np.zeros(6)
    cov[5][5] = 0.1 * 4 * np.pi**2 * g**3

    return a, adot, b, bdot, g, gdot, cov

def fit_observations(self):

    nobs = self.nobs
    ia = np.ones(6)

    if nobs < 2:
        raise ValueError('Not enough observations to fit an orbit for \
                          object {}'.format(self.hash))
    elif nobs == 2:
        fit_params = 4
    else:
        a, adot, b, bdot, g, gdot, cov = prelim_fit()
        if (cov[4][4] < 0) or cov[4][4]/g**2 > GDEGENERATE**2:
            fit_params = 4
        else: fit_params = 6

    if fit_params == 6:
