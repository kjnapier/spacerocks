class Ephemerides:

    '''
    Take in the ecliptic state vector (relative to the observer) and
    return any quantity
    '''

    def __init__(self, x, y, z, vx, vy, vz):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

        self.H = H


    @property
    def ra(self):
        return Angle(np.arctan2(self.y, self.x), u.rad).wrap_at(2 * np.pi * u.rad)

    @property
    def dec(self):
        return Angle(np.arcsin(self.z / sqrt(self.x**2 + self.y**2 + self.z**2)), u.rad)

    def ra_rate(self):
        return -(self.y * self.vx - self.x * self.vy) / (self.x**2 + self.y**2) * u.rad

    def dec_rate(self):
        return (-self.z * (self.x * self.vx + self.y * self.vy) + ((self.x**2 + self.y**2) * self.vz)) \
                / (sqrt(self.x**2 + self.y**2) * (self.x**2 + self.y**2 + self.z**2)) * u.rad

    def mag(self):
        pass

    def H(self):
        pass
