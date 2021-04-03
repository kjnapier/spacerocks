class Vector:

    def __init__(self, x, y, z):

        self.x = x
        self.y = y
        self.z = z

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

    @property
    def norm(self):
        return sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def unit(self):
        return self / self.norm 
