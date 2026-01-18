import numpy as np
from numpy.typing import NDArray
from manim import Polygon, Line3D


def poly_from_V(
    V: NDArray[np.float64],
    **style
) -> Polygon:
    return Polygon(*[tuple(p) for p in np.asarray(V)], **style)


def set_polygon_from_V(poly: Polygon, V: NDArray[np.float64]) -> Polygon:
    poly.set_points_as_corners([tuple(p) for p in np.asarray(V)])
    return poly


def line3d(p0, p1, **style) -> Line3D:
    return Line3D(p0, p1, **style)


def become_line3d(m: Line3D, p0, p1, **style) -> Line3D:
    m.become(Line3D(p0, p1, **style))
    return m
