import numpy as np
from typing import Tuple, List
from numpy.typing import NDArray


def make_square(
    side: float = 1.0,
    origin: Tuple[float, float] = (0.0, 0.0)
) -> Tuple[NDArray[np.float64], NDArray[np.int32]]:
    ox, oy = (float(origin[0]), float(origin[1]))
    s = float(side)
    V = np.array(
        [
            [ox + 0, oy + 0, 0.0],
            [ox + s, oy + 0, 0.0],
            [ox + s, oy + s, 0.0],
            [ox + 0, oy + s, 0.0],
        ],
        dtype=np.float64,
    )
    E = np.array([[0, 1], [1, 2], [2, 3], [3, 0]], dtype=np.int32)
    return V, E


def make_rect(
    width: float,
    height: float,
    origin: Tuple[float, float] = (0.0, 0.0)
) -> Tuple[NDArray[np.float64], NDArray[np.int32]]:
    ox, oy = (float(origin[0]), float(origin[1]))
    w = float(width)
    h = float(height)
    V = np.array(
        [
            [ox + 0, oy + 0, 0.0],
            [ox + w, oy + 0, 0.0],
            [ox + w, oy + h, 0.0],
            [ox + 0, oy + h, 0.0],
        ],
        dtype=np.float64,
    )
    E = np.array([[0, 1], [1, 2], [2, 3], [3, 0]], dtype=np.int32)
    return V, E


def make_equilateral_triangle(side: float):
    s = float(side)
    A = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    B = np.array([s, 0.0, 0.0], dtype=np.float64)
    C = np.array([0.5 * s, (np.sqrt(3) / 2) * s, 0.0], dtype=np.float64)
    V = np.array([A, B, C], dtype=np.float64)
    E = [(0, 1), (1, 2), (2, 0)]
    return V, E


def make_wall_from_edge(V_tri: NDArray[np.float64], edge, height: float):
    i, j = edge
    P0 = V_tri[i][:2]
    P1 = V_tri[j][:2]
    d = P1 - P0
    dx, dy = float(d[0]), float(d[1])
    length = float(np.hypot(dx, dy))
    if length == 0.0:
        raise ValueError("Degenerate edge for wall")

    # outward unit normal = right of directed edge
    u_out = np.array([dy, -dx], dtype=np.float64) / length
    Q0 = P0 + u_out * float(height)
    Q1 = P1 + u_out * float(height)

    P0_3 = np.array([P0[0], P0[1], 0.0], dtype=np.float64)
    P1_3 = np.array([P1[0], P1[1], 0.0], dtype=np.float64)
    Q1_3 = np.array([Q1[0], Q1[1], 0.0], dtype=np.float64)
    Q0_3 = np.array([Q0[0], Q0[1], 0.0], dtype=np.float64)

    V_wall = np.array([P0_3, P1_3, Q1_3, Q0_3], dtype=np.float64)
    E_wall = [(0, 1), (1, 2), (2, 3), (3, 0)]
    return V_wall, E_wall


def make_lid_triangle_from_wall(V_wall0: NDArray[np.float64], side: float):
    # Wall vertices: [P0, P1, Q1, Q0], outer edge = (2,3) = Q1 -> Q0
    Q1 = V_wall0[2][:2]
    Q0 = V_wall0[3][:2]

    L0 = Q0
    L1 = Q1
    d = L1 - L0
    dx, dy = float(d[0]), float(d[1])
    length = float(np.hypot(dx, dy))
    if length == 0.0:
        raise ValueError("Degenerate lid base edge")

    side_eff = length  # trust actual geometry

    # choose interior on "outside" of the wall
    u_int = -np.array([-dy, dx], dtype=np.float64) / length
    h_tri = (np.sqrt(3) / 2) * side_eff
    L2 = (L0 + L1) / 2 + u_int * h_tri

    L0_3 = np.array([L0[0], L0[1], 0.0], dtype=np.float64)
    L1_3 = np.array([L1[0], L1[1], 0.0], dtype=np.float64)
    L2_3 = np.array([L2[0], L2[1], 0.0], dtype=np.float64)

    V_lid = np.array([L0_3, L1_3, L2_3], dtype=np.float64)
    E_lid = [(0, 1), (1, 2), (2, 0)]
    lid_hinge = (0, 1)
    return V_lid, E_lid, lid_hinge
