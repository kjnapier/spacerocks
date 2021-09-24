'''
Do a Gauss orbit fit.
'''
from numpy import sin, cos, roots, isreal, array

from .vector import Vector
from .spacerock import SpaceRock
from .constants import mu_bary as mu
from .constants import epsilon
mu = mu.value


def gauss(ras, decs, epochs, observer) -> SpaceRock:

    # Step -1: Calculate the xyz position of the observer
    R1 = Vector(observer[0].x.au, observer[0].y.au, observer[0].z.au)
    R2 = Vector(observer[1].x.au, observer[1].y.au, observer[1].z.au)
    R3 = Vector(observer[2].x.au, observer[2].y.au, observer[2].z.au)

    # Step 0: Calculate unit vectors
    xhat = cos(decs) * cos(ras)
    yhat = cos(decs) * sin(ras)
    zhat = sin(decs)

    rho1 = Vector(xhat[0], yhat[0], zhat[0])
    rho2 = Vector(xhat[1], yhat[1], zhat[1])
    rho3 = Vector(xhat[2], yhat[2], zhat[2])

    # Step 1: calculate time intervals
    t1, t2, t3 = epochs.tdb.jd
    tau1 = t1 - t2
    tau3 = t3 - t2
    tau = t3 - t1

    # Step 2: Calculate cross products, take the cross products of the observational unit direction
    p1 = rho2.cross(rho3)
    p2 = rho1.cross(rho3)
    p3 = rho1.cross(rho2)

    # Step 3: Compute a common scalar quantity: The triple product
    D0 = rho1.dot(p1)

    # Step 4: Compute the D tensor
    D = array([[R1.dot(p1), R1.dot(p2), R1.dot(p3)],
               [R2.dot(p1), R2.dot(p2), R2.dot(p3)],
               [R3.dot(p1), R3.dot(p2), R3.dot(p3)]])

    # Step 5: Calculate scalar position coefficients
    A = (1/D0) * (-D[0][1] * (tau3/tau) + D[1][1] + D[2][1] * (tau1/tau))
    B = (1/(6 * D0)) * (D[0][1] * (tau3**2 - tau**2) * (tau3/tau) + D[2][1] * (tau**2 - tau1**2) * (tau1/tau))
    E = R2.dot(rho2)

    # Step 6: Calculate the squared scalar distance of the second observation
    R2sq = R2.dot(R2)

    # Step 7: Calculate the coefficients of the scalar distance polynomial for the second observation of the orbiting body
    a = -(A**2 + 2 * A * E + R2sq)
    b = -2 * mu * B * (A + E)
    c = -mu**2 * B**2

    # Step 8: Find the root of the scalar distance polynomial for the second observation of the orbiting body
    coeff = [1, 0, a.value, 0, 0, b.value, 0, 0, c.value]
    rts = roots(coeff)
    rts = rts[isreal(rts)]
    rts = rts[rts > 0]

    x_out = []
    y_out = []
    z_out = []
    vx_out = []
    vy_out = []
    vz_out = []
    epoch_out = []

    for root in rts:

        root = root.real

        # Step 9: Calculate the slant range, the distance from the observer point to the orbiting body at their respective time
        a1 = (1/D0) * ((6 * (D[2][0] * (tau1/tau3) + D[1][0] * (tau/tau3)) * root**3 + mu * D[2][0] * (tau**2 - tau1**2) * (tau1/tau3)) / (6 * root**3 + mu * (tau**2 - tau3**2)) - D[0][0])
        a2 = A + (mu * B) / root**3
        a3 = (1/D0) * ((6 * (D[0][2] * (tau3/tau1) - D[1][2] * (tau/tau1)) * root**3 + mu * D[0][2] * (tau**2 - tau3**2) * (tau3/tau1)) / (6 * root**3 + mu * (tau**2 - tau1**2)) - D[2][2])

        # Step 10: Calculate the orbiting body position vectors, by adding the observer position vector to the slant direction vector


        r1 = array([R1.x, R1.y, R1.z]) + a1 * array([rho1.x, rho1.y, rho1.z])
        r2 = array([R2.x, R2.y, R2.z]) + a2 * array([rho2.x, rho2.y, rho2.z])
        r3 = array([R3.x, R3.y, R3.z]) + a3 * array([rho3.x, rho3.y, rho3.z])

        # Step 11: Calculate the Lagrange Coefficients
        f1 = 1 - 0.5 * (mu/root**3) * tau1**2
        f3 = 1 - 0.5 * (mu/root**3) * tau3**2
        g1 = tau1 - (1/6) * (mu / root**3) * tau1**3
        g3 = tau3 - (1/6) * (mu / root**3) * tau3**3

        # Step 12: Calculate the velocity vector for the second observation of the orbiting body
        v2 = (-f3 * r1 + f1 * r3) / (f1 * g3 - f3 * g1)

        # Step 13: The orbital state vectors have now been found, the position (r2) and velocity (v2) vector for the second observation of the orbiting body. With these two vectors, the orbital elements can be found and the orbit determined
        x, y, z = r2
        vx, vy, vz = v2

        y_ec = y * cos(epsilon) + z * sin(epsilon)
        z_ec = -y * sin(epsilon) + z * cos(epsilon)
        vy_ec = vy * cos(epsilon) + vz * sin(epsilon)
        vz_ec = -vy * sin(epsilon) + vz * cos(epsilon)

        x_out.append(x.value)
        y_out.append(y_ec.value)
        z_out.append(z_ec.value)
        vx_out.append(vx.value)
        vy_out.append(vy_ec.value)
        vz_out.append(vz_ec.value)
        epoch_out.append(t2)

    return x_out, y_out, z_out, vx_out, vy_out, vz_out, epoch_out