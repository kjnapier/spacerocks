#import emcee

from .spacerock import SpaceRock
from .gauss import gauss
from .observer import Observer


class OrbitFitter:

    def __init__(self, ra, dec, epoch, obscode):

        self.ra = ra
        self.dec = dec
        self.epoch = epoch
        self.obscode = obscode

        #self.observer = observer(obscode=self.obscode, epoch=self.epoch)
        self.observer = Observer(self.obscode, self.epoch)

    def gauss_guess(self):
        x, y, z, vx, vy, vz, epoch = gauss(self.ra, self.dec, self.epoch, self.observer)
        return SpaceRock(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=epoch)    

    # def mcmc(self):

    #     def log_prior(theta):
    #         a, e, inc, arg, node, true_anomaly = theta
    #         if a > 0 and 0 <= e < 1 and 0.0 < inc < 180.0 and 0.0 < arg < 360 and 0.0 < node < 360 and 0.0 < true_anomaly < 360:
    #             return 0.0
    #         elif a < 0 and e > 1 and 0.0 < inc < 180.0 and 0.0 < arg < 360 and 0.0 < node < 360 and 0.0 < true_anomaly < 360:
    #             return 0.0
    #         else:
    #             return -np.inf

    #     def log_likelihood(theta, obs_ra, obs_dec, t, obscode):

    #         a, e, inc, arg, node, true_anomaly = theta
    #         ucty = 0.15/3600
    #         rock = SpaceRock(a=a, e=e, inc=inc, arg=arg, node=node, true_anomaly=true_anomaly, epoch=t[0])
    #         prop, _, _ = rock.propagate(epochs=t, model=1)
    #         #prop = rock.analytic_propagate(epochs=t, propagate_frame='barycentric')
    #         obs = prop.observe(obscode=obscode)

    #         ra = obs.ra.deg
    #         dec = obs.dec.deg
    #         ra_rate = obs.ra_rate.value
    #         dec_rate = obs.dec_rate.value

    #         csq = np.sum((ra - obs_ra)**2 / ucty**2 + (dec - obs_dec)**2 / ucty**2)

    #     def log_probability(theta, obs_ra, obs_dec, t, obscode):
    #         lp = log_prior(theta)
    #         if not np.isfinite(lp):
    #             return -np.inf
    #         return lp + log_likelihood(theta, obs_ra, obs_dec, t, obscode)
    #         pass

    #     args = [self.ra.deg, self.dec.deg, self.epoch.tdb.jd, self.obscode]
    #     sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=args)

    #     start = self.gauss_guess()
    #     p = np.array([start.a.au,
    #                   start.e,
    #                   start.inc.deg,
    #                   start.arg.deg,
    #                   start.node.deg,
    #                   start.true_anomaly.deg]).T

    #     pos, prob, state = sampler.run_mcmc(p, 100, skip_initial_state_check=True, progress=True)
