import numpy as np


def flat_to_cylinder_coord(p_flat, R: float, L: float):
    x, y, z = p_flat
    u = x / (L / 2)              # [-1,1]
    theta = np.pi * (1.0 + u)    # [-1->0], [0->pi], [1->2pi]
    Xc = R * np.cos(theta)
    Zc = R * np.sin(theta)
    return np.array([Xc, y, Zc], dtype=np.float64)


def map_flat_to_cylinder(p_flat, t: float, R: float, L: float):
    p_flat = np.asarray(p_flat, dtype=np.float64)
    p_cyl = flat_to_cylinder_coord(p_flat, R=R, L=L)
    return (1 - t) * p_flat + t * p_cyl
