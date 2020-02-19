# Take a date (in Universal Time), and return the it in the desired form.
# Vectorize this!
import numpy as np
class Date:

    def __init__(self, day, month, year, UT):

        self.day = day
        self.month = month
        self.year = year
        self.UT = UT # Universal Time
        self.JD = np.vectorize(self.jd)
        
    def jd(self):

        if self.month > 2:
            y = self.year
            m = self.month
        else:
            y = self.year - 1
            m = self.month + 12

        if self.year > 1582:
            B = int(y / 400) - int(y / 100)
        elif self.year < 1582:
            B = -2
        else:
            if self.month < 10:
                B = -2
            elif self.month > 10:
                B = int(y / 400) - int(y /100)
            else:
                if self.day <= 4:
                    B = -2
                elif self.day >= 15:
                    B = int(y / 400) - int(y / 100)

        return int(365.25 * y) + int(30.6001 * (m + 1)) + 1720996.5 + B + self.day + self.UT / 24

    def MJD(self):
        return self.JD() - 2400000.5

print(Date(np.array([1, 1]), np.array([1, 1]), np.array([1, 1]), np.array([1, 1])).JD())