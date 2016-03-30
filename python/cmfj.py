"""
Climate model from Advanced Atmospheric Physics lecture series, written in python
"""
#pylint: disable = invalid-name

import numpy as np
#from numba import jit, vectorize, float64
import datetime

CONST_G = 9.81
CONST_CP = 1.006e3


def schwarzschild(dtau, w0, g, B, Ag, Eup, Edn, divE):
    """
    Compute radiative transfer according to schwarzschild equation, neglecting scattering
    returns Irradiace downward, upward and divergence of radiation
    """
    Nmu = 5
    dmu = 1./Nmu

    ke = len(dtau)

    kabs = dtau * (1. - w0)

    Eup[:] = 0
    Edn[:] = 0

    for imu in xrange(Nmu):
        mu = (imu + .5)*dmu

        T = np.exp(-kabs / mu)

        Lup = B[ke-1]*(1.-Ag)
        Eup[ke] = Eup[ke] + Lup*mu

        Ldn = 0
        Edn[0] = Edn[0] + Ldn*mu

        for k in xrange(ke-1, -1, -1):
            Lup = Lup * T[k] + B[k] * (1.-T[k])
            Eup[k] = Eup[k] + Lup*mu

        for k in xrange(0, ke):
            Ldn = Ldn * T[k] + B[k] * (1.-T[k])
            Edn[k+1] = Edn[k+1] + Ldn*mu

    Eup[:] = Eup[:] * 2 * np.pi * dmu
    Edn[:] = Edn[:] * 2 * np.pi * dmu
    divE[:] = Eup[1:] - Eup[:-1] + Edn[:-1] - Edn[1:]
    divE[-1] += Edn[-1] - Eup[-1]

#jitschwarzschild = jit('void( float64[:], float64[:], float64[:], float64[:], float64, float64[:], float64[:], float64[:])', nopython=True, nogil=True)(schwarzschild)

#@vectorize([float64(float64)])
def stephanboltzmann(T):
    """ return blackbody radiation according to Stephan Boltzmanns law """
    return 5.67e-8 * T**4 /np.pi

#@vectorize([float64(float64, float64, float64)])
def T2theta(T, p, p0):
    """ convert Temperature in K to potential Temperature at reference pressure p0 """
    return T * (p0/p)**(2./7)

#@vectorize([float64(float64, float64, float64)])
def theta2T(T, p, p0):
    """ convert potential Temperature to Temperature """
    return T / (p0/p)**(2./7)

#@jit('float64(float64[:], float64[:])', nopython=True)
def rmse(a, b):
    return np.sqrt(np.mean((a-b)**2))

#@profile
def cmfj(nlyr=100, T0=288, p0=1013e2):
    """
    Main function for climate model run
    """
    a = datetime.datetime.now()

    pt = np.linspace(0, p0, nlyr+1)
    pm = (pt[1:] + pt[:-1])/2.
    T = (-np.arange(nlyr)[::-1] + T0).astype(np.float64)
    lastT = T.copy()
    divE = np.zeros(nlyr)

    converged = False
    while not converged:
        dtau = 1. * np.ones(nlyr) / nlyr
        w0 = np.zeros(nlyr)
        g = np.zeros(nlyr)
        B = stephanboltzmann(T)

        Eup = np.zeros(nlyr+1)
        Edn = np.zeros(nlyr+1)
        divE = np.zeros(nlyr)

        schwarzschild(dtau, w0, g, B, 0, Eup, Edn, divE)
#        jitschwarzschild(dtau, w0, g, B, np.float64(0), Eup, Edn, divE)

        divE[-1] += 235

        Tinc = divE *CONST_G / CONST_CP / (pt[1:]-pt[:-1])

        dt = 1. / np.max(np.abs(Tinc))

        T = T + Tinc*dt

        new_theta = np.sort(T2theta(T, pm, p0))[::-1]
        T = theta2T(new_theta, pm, p0)

        residual = rmse(lastT, T)
        print residual
        if residual < 1e-4:
            converged = True
        lastT = T.copy()

    b = datetime.datetime.now()

    delta = b - a

    print "== Converged =="
    for k in xrange(nlyr):
        print k, " Pressure %s  Temperature %s" % (pm[k], T[k])

    print "== Time: "+ str(delta.total_seconds() *1000)+ " ms"

if __name__ == "__main__":
    cmfj()
