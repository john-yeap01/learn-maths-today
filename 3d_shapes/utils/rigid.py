import numpy as np
from typing import Iterable, Tuple
from numpy.typing import NDArray


def unit_vec(v: Iterable[float]) -> NDArray[np.float64]:
    x = np.asarray(v, dtype=float).reshape(-1)
    if x.size != 3:
        raise ValueError("Expected 3 components")
    n = float(np.linalg.norm(x))
    if n == 0.0:
        raise ValueError("zero length vector")
    return (x / n).astype(np.float64, copy=False)


def skew(v: Iterable[float]) -> NDArray[np.float64]:
    vec = np.asarray(v, dtype=float).reshape(-1)
    if vec.size != 3:
        raise ValueError("Expected 3 components")
    x_, y_, z_ = (float(vec[0]), float(vec[1]), float(vec[2]))
    return np.array(
        [[0, -z_, y_],
         [z_, 0, -x_],
         [-y_, x_, 0]],
        dtype=np.float64,
    )


def rodrigues(axis: Iterable[float], theta: float) -> NDArray[np.float64]:
    k = unit_vec(axis)
    K = skew(k)
    s = float(np.sin(theta))
    c = float(np.cos(theta))
    I = np.eye(3, dtype=np.float64)
    return I + s * K + (1 - c) * (K @ K)


def se3_about_edge(
    p0: NDArray[np.float64],
    p1: NDArray[np.float64],
    theta: float
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    p0 = np.asarray(p0, dtype=np.float64).reshape(3)
    p1 = np.asarray(p1, dtype=np.float64).reshape(3)
    R = rodrigues(p1 - p0, theta)
    t = p0 - (R @ p0)  # keep the whole hinge line fixed
    return R, t


def apply_se3(
    V: NDArray[np.float64],
    R: NDArray[np.float64],
    t: NDArray[np.float64]
) -> NDArray[np.float64]:
    V = np.asarray(V, dtype=np.float64)
    return (R @ V.T).T + t


def compose_se3(
    R1: NDArray[np.float64],
    t1: NDArray[np.float64],
    R2: NDArray[np.float64],
    t2: NDArray[np.float64]
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    # x' = R2(R1 x + t1) + t2 = (R2 R1) x + (R2 t1 + t2)
    R = R2 @ R1
    t = (R2 @ t1) + t2
    return R, t
