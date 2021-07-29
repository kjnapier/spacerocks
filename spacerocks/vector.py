from numpy import sqrt, sin, cos
import copy

class Vector:

    def __init__(self, x, y, z):

        self.x = x
        self.y = y
        self.z = z

    def __getitem__(self, idx):
        '''
        This method allows you to index a Vector object.
        '''
        p = copy.copy(self)
        for attr in self.__dict__.keys():
            setattr(p, attr, getattr(self, attr)[idx])

        return p

    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __rmul__(self, c):
        return Vector(c * self.x, c * self.y, c * self.z)

    def __truediv__(self, c):
        return Vector(self.x / c, self.y / c, self.z / c)

    def dot(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other):
        x = self.y * other.z - self.z * other.y
        y = self.z * other.x - self.x * other.z
        z = self.x * other.y - self.y * other.x
        return Vector(x, y, z)

    def euler_rotation(self, a, b, c):
        sa = sin(a)
        sb = sin(b)
        sc = sin(c)
        ca = cos(a)
        cb = cos(b)
        cc = cos(c)
        # Currently assumes z=0. I'm going to modify this to the genreal case.
        xrot = self.x * (ca * cc - sa * sc * cb) - self.y * (sa * cc + ca * sc * cb)
        yrot = self.x * (ca * sc + sa * cc * cb) + self.y * (ca * cc * cb - sa * sc)
        zrot = self.x * (sa * sb) + self.y * (ca * sb)
        #xrot = self.x * (cos(a) * cos(c) - sin(a) * sin(c) * cos(b)) - self.y * (sin(a) * cos(c) + cos(a) * sin(c) * cos(b))
        #yrot = self.x * (cos(a) * sin(c) + sin(a) * cos(c) * cos(b)) + self.y * (cos(a) * cos(c) * cos(b) - sin(a) * sin(c))
        #zrot = self.x * (sin(a) * sin(b)) + self.y * (cos(a) * sin(b))
        return Vector(xrot, yrot, zrot)


    @property
    def norm(self):
        return sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def unit(self):
        return self / self.norm
