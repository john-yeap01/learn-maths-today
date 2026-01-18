import numpy as np
from numpy.typing import NDArray
from manim import Polygon, Line3D, Dot3D, ValueTracker
from .rigid import se3_about_edge, apply_se3
from .manim_bridge import set_polygon_from_V, become_line3d


def make_face_hinge_updater(
    V0: NDArray[np.float64],
    p0: NDArray[np.float64],
    p1: NDArray[np.float64],
    tracker: ValueTracker,
    sign: float = +1.0,
):
    V0 = np.asarray(V0, dtype=np.float64).copy()
    p0 = np.asarray(p0, dtype=np.float64).reshape(3)
    p1 = np.asarray(p1, dtype=np.float64).reshape(3)

    def _upd(poly: Polygon):
        R, t = se3_about_edge(p0, p1, float(sign) * tracker.get_value())
        V = apply_se3(V0, R, t)
        return set_polygon_from_V(poly, V)

    return _upd


def make_point_hinge_updater(
    x0: NDArray[np.float64],
    p0: NDArray[np.float64],
    p1: NDArray[np.float64],
    tracker: ValueTracker,
    sign: float = +1.0,
):
    x0 = np.asarray(x0, dtype=np.float64).reshape(3)
    p0 = np.asarray(p0, dtype=np.float64).reshape(3)
    p1 = np.asarray(p1, dtype=np.float64).reshape(3)

    def _upd(dot: Dot3D):
        R, t = se3_about_edge(p0, p1, float(sign) * tracker.get_value())
        dot.move_to(R @ x0 + t)
        return dot

    return _upd


def make_child_hinge_updater(
    V_child0: NDArray[np.float64],
    parent_p0: NDArray[np.float64],
    parent_p1: NDArray[np.float64],
    parent_tracker: ValueTracker,
    parent_sign: float,
    child_hinge_idx: tuple[int, int],
    child_tracker: ValueTracker,
    child_sign: float,
):
    V_child0 = np.asarray(V_child0, dtype=np.float64).copy()
    parent_p0 = np.asarray(parent_p0, dtype=np.float64).reshape(3)
    parent_p1 = np.asarray(parent_p1, dtype=np.float64).reshape(3)
    h0, h1 = int(child_hinge_idx[0]), int(child_hinge_idx[1])

    def _upd(poly: Polygon):
        # parent
        R_p, t_p = se3_about_edge(parent_p0, parent_p1, float(parent_sign) * parent_tracker.get_value())
        Vp = apply_se3(V_child0, R_p, t_p)

        # moved hinge (world)
        H0, H1 = Vp[h0], Vp[h1]

        # child
        R_c, t_c = se3_about_edge(H0, H1, float(child_sign) * child_tracker.get_value())
        V = apply_se3(Vp, R_c, t_c)

        return set_polygon_from_V(poly, V)

    return _upd


def make_line_child_hinge_updater(
    V_child0: NDArray[np.float64],
    parent_p0: NDArray[np.float64],
    parent_p1: NDArray[np.float64],
    parent_tracker: ValueTracker,
    parent_sign: float,
    line_idx: tuple[int, int],
    child_tracker: ValueTracker | None = None,
    child_sign: float = +1.0,
):
    """
    Updates a Line3D that should track the child hinge.
    If child_tracker is None => just apply parent transform (no second rotation).
    If child_tracker provided => apply parent then child around the moved hinge.
    """
    V_child0 = np.asarray(V_child0, dtype=np.float64).copy()
    parent_p0 = np.asarray(parent_p0, dtype=np.float64).reshape(3)
    parent_p1 = np.asarray(parent_p1, dtype=np.float64).reshape(3)
    i0, i1 = int(line_idx[0]), int(line_idx[1])

    def _upd(line: Line3D):
        R_p, t_p = se3_about_edge(parent_p0, parent_p1, float(parent_sign) * parent_tracker.get_value())
        Vp = apply_se3(V_child0, R_p, t_p)
        H0, H1 = Vp[i0], Vp[i1]

        if child_tracker is not None:
            R_c, t_c = se3_about_edge(H0, H1, float(child_sign) * child_tracker.get_value())
            H0 = R_c @ H0 + t_c
            H1 = R_c @ H1 + t_c

        return become_line3d(line, H0, H1, color=line.color)

    return _upd
