import numpy as np
import angles

precise = True

class SpaceRock:
    
    def __init__(self, a, e, i, node, omega, tau, epoch):

        self.a = a
        self.e = e
        self.i = i
        self.node = node
        self.omega = omega
        self.epoch = epoch
        self.tau = tau # obs_date
        self.T = self.a**(3/2)
        self.n = (2*np.pi)/self.T
        self.M = self.n * (self.epoch - self.tau)
        
    if precise == True:

        def E(self):
            calc_E = self.M
            for kk in range(100):
                calc_E = self.M + self.e * np.sin(calc_E)
            return calc_E

    else:

        def cal_E(self, e, M):
            # compute eccentric anomaly E
            f = lambda calc_E, M, e: calc_E - e * np.sin(E) - M
            E0 = M
            calc_E = newton(f, E0, args=(M, e))
            return calc_E


    
    def x(self):
        return self.a * (np.cos(self.E()) - self.e)
    
    def y(self):
        return self.a * np.sqrt(1-self.e**2) * np.sin(self.E())
    
    def z(self):
        return 0;
    
    def xyz(self):
        return np.array([[self.x()],[self.y()],[self.z()]])
    
    def R1(self):
        return np.array([[np.cos(self.omega), -np.sin(self.omega), 0],
                         [np.sin(self.omega), np.cos(self.omega), 0],
                         [0, 0, 1]]);
    
    def R2(self):
        return np.array([[1, 0, 0],
                         [0, np.cos(self.i), -np.sin(self.i)],
                         [0, np.sin(self.i), np.cos(self.i)]]);
    
    def R3(self):
        return np.array([[np.cos(self.node), -np.sin(self.node), 0],
                         [np.sin(self.node), np.cos(self.node), 0],
                         [0, 0, 1]]);
    
    def RMatrix(self):
        return multi_dot([self.R3(), self.R2(), self.R1()])
    
    def XYZ(self):
        return np.dot(self.RMatrix(), self.xyz())
    
    def R(self):
        return np.sqrt(np.sum(self.XYZ() * self.XYZ()))

class Date:

    def __init__(self, day, month, year, UT):
    
        self.day = day
        self.month = month
        self.year = year
        self.UT = UT # Universal Time
        self.JD = np.vectorize(self.jd)
    
    def test(self):
        for ii in range(len(self.day)):
            if self.day[ii] > 1:
                return self.day[ii] * 10
    
    def jd(self):
    
        if self.month > 2:
            y = self.year
            m = self.month
        else:
            y = self.year - 1
            m = self.month + 12
    
        if self.year > 1582:
            B = y // 400 - y // 100
        elif self.year < 1582:
            B = -2
        else:
            if self.month < 10:
                B = -2
            elif self.month > 10:
                B = y // 400 - y // 100
            else:
                if self.day <= 4:
                    B = -2
                elif self.day >= 15:
                    B = y // 400 - y // 100
    
        return np.array(int(365.25 * y) + int(30.6001 * (m + 1)) + 1720996.5 + B + self.day + self.UT / 24)
    
    def MJD(self):
        return self.JD() - 2400000.5
    
