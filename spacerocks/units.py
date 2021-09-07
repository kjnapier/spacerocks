from astropy import units as u

# Units manager for spacerocks.
# Usage:
# from astropy import units as u
# units = Units()
# units.angle = u.rad


class Units:

    def __init__(self):
        self.set_default()

    def set_default(self):
        self.distance = u.au
        self.angle = u.deg
        
        self.timescale = 'utc'
        self.timeformat = None
        self.speed = u.au / u.d

        self.size = u.km
        self.density = u.g / u.cm**3
        self.mass = u.kg

        self.ra = u.deg
        self.dec = u.deg

    def current(self):
        print("{:<20} {:<15}".format('Quantity', 'Unit'))
        print('---------------------------------------')
        for k, v in self.__dict__.items():
            print("{:<20} {:<15}".format(k, str(v)))
