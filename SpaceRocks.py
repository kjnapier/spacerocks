import numpy as numpy
import angles
import dates

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