"""STEP 10 — LV models linking end-diastolic geometry, strain and EF (lengths in mm, volumes in mm^3, strains as fractions,
negative = shortening). Model estimates only: the meshes are end-diastolic, so no strain is measured.

Our model (STEP 10.1-10.3): half prolate ellipsoid, uniform wall, strain at midwall, incompressible myocardium.
Stokke model (STEP 10.6): Equation 1 of Stokke et al., JACC 2017;70:942-954 (doi:10.1016/j.jacc.2017.06.046), transcribed
from the article: strain at a layer f x wall thickness out from the endocardium (default f = 0.33), myocardial volume
shrinking to c x its diastolic value in systole (default c = 0.95); the long axis cancels out.
"""
import numpy as np
from scipy import optimize

K = 2 / 3 * np.pi  # half-ellipsoid volume = K * a^2 * c


# STEP 10.1 — geometry from EDV, myocardial volume and LV length
def geometry(edv, myo, L):
    r = np.sqrt(edv / (K * L))
    T = optimize.brentq(lambda t: K * ((r + t) ** 2 * (L + t) - r * r * L) - myo, 0, 50)
    return r, T


# STEP 10.2 — EF for given strains
def ef(r, L, T, gls, gcs):
    myo = K * ((r + T) ** 2 * (L + T) - r * r * L)
    rm, lm = (r + T / 2) * (1 + gcs), (L + T / 2) * (1 + gls)
    f = lambda t: K * ((rm + t / 2) ** 2 * (lm + t / 2) - (rm - t / 2) ** 2 * (lm - t / 2)) - myo  # noqa: E731
    if f(2 * min(rm, lm)) < 0:  # wall cannot hold the myocardium: cavity obliterated
        return 1.0
    t = optimize.brentq(f, 0, 2 * min(rm, lm))
    return 1 - (rm - t / 2) ** 2 * (lm - t / 2) / (r * r * L)


# STEP 10.3 — GLS needed for a target EF, with GCS = ratio x GLS
def gls_needed(r, L, T, target, ratio):
    return optimize.brentq(lambda s: ef(r, L, T, s, ratio * s) - target, -0.6, 0)


# STEP 10.6 — Stokke Equation 1: d = ED short-axis (cavity) diameter, w = ED wall thickness
def stokke_ef(d, w, gls, gcs, f=0.33, c=0.95):
    R, fw = d / 2, f * w
    return 1 - ((1 + gcs) ** 2 * (1 + gls) * (R + fw) ** 2 - 2 * c * fw * (R + fw) + c * fw * fw) / R ** 2


def stokke_gls_needed(d, w, target, ratio, f=0.33, c=0.95):
    return optimize.brentq(lambda s: stokke_ef(d, w, s, ratio * s, f, c) - target, -0.6, 0)
