import numpy as np

from .detection import Detection
from ..units import Units
from ..spacerock import SpaceRock
from ..constants import mu_bary as mu
mu = mu.value

def gauss(detections: list[Detection], units: Units = Units()) -> SpaceRock:

    # observer = observer.change_frame('J2000')

    # # Step -1: Calculate the xyz position of the observer
    # R1 = Vector(observer.x[0].au, observer.y[0].au, observer.z[0].au)
    # R2 = Vector(observer.x[1].au, observer.y[1].au, observer.z[1].au)
    # R3 = Vector(observer.x[2].au, observer.y[2].au, observer.z[2].au)

    # # Step 0: Calculate unit vectors
    # rho1 = detections[0].pointing_vector
    # rho2 = detections[1].pointing_vector
    # rho3 = detections[2].pointing_vector

    R1, R2, R3 = [detection.observer.position for detection in detections]
    rho1, rho2, rho3 = [detection.pointing_vector for detection in detections]

    # Step 1: calculate time intervals
    t1, t2, t3 = (detection.epoch for detection in detections)
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
    D = np.array([[R1.dot(p1), R1.dot(p2), R1.dot(p3)],
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
    coeff = [1, 0, a, 0, 0, b, 0, 0, c]

    try:
        rts = np.roots(coeff)
    except Exception as e:
        return None

    rts = rts[np.isreal(rts)]
    rts = rts[rts > 0]
    rts = [root.real for root in rts]

    x_out = []
    y_out = []
    z_out = []
    vx_out = []
    vy_out = []
    vz_out = []
    epoch_out = []

    for root in rts:

        # Step 9: Calculate the slant range, the distance from the observer point to the orbiting body at their respective time
        a1 = (1/D0) * ((6 * (D[2][0] * (tau1/tau3) + D[1][0] * (tau/tau3)) * root**3 + mu * D[2][0] * (tau**2 - tau1**2) * (tau1/tau3)) / (6 * root**3 + mu * (tau**2 - tau3**2)) - D[0][0])
        a2 = A + (mu * B) / root**3
        a3 = (1/D0) * ((6 * (D[0][2] * (tau3/tau1) - D[1][2] * (tau/tau1)) * root**3 + mu * D[0][2] * (tau**2 - tau3**2) * (tau3/tau1)) / (6 * root**3 + mu * (tau**2 - tau1**2)) - D[2][2])

        # Step 10: Calculate the orbiting body position vectors, by adding the observer position vector to the slant direction vector
        r1 = np.array([R1.x, R1.y, R1.z]) + a1 * np.array([rho1.x, rho1.y, rho1.z])
        r2 = np.array([R2.x, R2.y, R2.z]) + a2 * np.array([rho2.x, rho2.y, rho2.z])
        r3 = np.array([R3.x, R3.y, R3.z]) + a3 * np.array([rho3.x, rho3.y, rho3.z])

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

        x_out.append(x)
        y_out.append(y)
        z_out.append(z)
        vx_out.append(vx)
        vy_out.append(vy)
        vz_out.append(vz)
        epoch_out.append(t2)

    rock = SpaceRock(x=x_out, y=y_out, z=z_out, vx=vx_out, vy=vy_out, vz=vz_out, epoch=epoch_out, origin='ssb', frame='J2000', units=units)

    return rock