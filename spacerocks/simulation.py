import rebound
import os
import spiceypy as spice
import pkg_resources
from astropy.time import Time

from .spice import SpiceBody

SPICE_PATH = pkg_resources.resource_filename('spacerocks', 'data/spice')

class Simulation(rebound.Simulation):

    def __init__(self, **kwargs):
        pass

    def set_simulation(self, startdate, model):

        sun = SpiceBody(spiceid='Sun')
        mercury = SpiceBody(spiceid='Mercury Barycenter')
        venus = SpiceBody(spiceid='Venus Barycenter')
        earth = SpiceBody(spiceid='Earth')
        moon = SpiceBody(spiceid='Moon')
        mars = SpiceBody(spiceid='Mars Barycenter')
        jupiter = SpiceBody(spiceid='Jupiter Barycenter')
        saturn = SpiceBody(spiceid='Saturn Barycenter')
        uranus = SpiceBody(spiceid='Uranus Barycenter')
        neptune = SpiceBody(spiceid='Neptune Barycenter')
        pluto = SpiceBody(spiceid='Pluto Barycenter')

        if model == 0:
            active_bodies = [sun]
            names = ['Sun']
    
        elif model == 1:
            active_bodies = [sun, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']

        elif model == 2:
            active_bodies = [sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune, pluto]
            names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']

        elif model == 3:

            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000001.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000002.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000003.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000004.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000007.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000010.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000015.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000016.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000031.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000052.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000065.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000087.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000088.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000107.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000511.bsp'))
            spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000704.bsp'))

            ceres = SpiceBody(spiceid='Ceres')
            vesta = SpiceBody(spiceid='Vesta')
            pallas = SpiceBody(spiceid='Pallas')
            juno = SpiceBody(spiceid='2000003')
            iris = SpiceBody(spiceid='2000007')
            hygiea = SpiceBody(spiceid='2000010')
            eunomia = SpiceBody(spiceid='2000015')
            psyche = SpiceBody(spiceid='2000016')
            euphrosyne = SpiceBody(spiceid='2000031')
            europa = SpiceBody(spiceid='2000052')
            cybele = SpiceBody(spiceid='2000065')
            sylvia = SpiceBody(spiceid='2000087')
            thisbe = SpiceBody(spiceid='2000088')
            camilla = SpiceBody(spiceid='2000107')
            davida = SpiceBody(spiceid='2000511')
            interamnia = SpiceBody(spiceid='2000704')
            
            active_bodies = [sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune, pluto, 
                             ceres, vesta, pallas, interamnia, juno, camilla, iris, hygiea, eunomia, psyche, 
                             euphrosyne, europa, cybele, sylvia, thisbe, davida]
            names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 
                     'Ceres', 'Vesta', 'Pallas', 'Interamnia', 'Juno', 'Camilla', 'Iris', 'Hygiea', 'Eunomia', 
                     'Psyche', 'Euphrosyne', 'Europa', 'Cybele', 'Sylvia', 'Thisbe', 'Davida']

        else:
            raise ValueError('Model not recognized. Check the documentation.')

        startdate = Time(startdate, scale='tdb', format='jd')

        bodies = [body.at(startdate) for body in active_bodies]

        sim = rebound.Simulation()
        sim.units = ('day', 'AU', 'Msun')

        for body, name, ab in zip(bodies, names, active_bodies):
            sim.add(x=body.x.au, 
                    y=body.y.au, 
                    z=body.z.au,
                    vx=body.vx.au, 
                    vy=body.vy.au, 
                    vz=body.vz.au,
                    m=ab.mass.value, 
                    hash=name)

        sim.N_active = len(bodies)

        #if gr == True:
        #    rebx = reboundx.Extras(sim)
        #    gr = rebx.load_force('gr_full')
        #    gr.params["c"] = constants.C
        #    rebx.add_force(gr)

        sim.testparticle_type = 0
        sim.integrator = 'ias15'
        #sim.integrator = 'mercurius'
        

        return sim, names