import numpy as np
from numpy.typing import NDArray
from manim import VGroup, Text, Line3D


def vertex_labels(V: NDArray[np.float64], font_size=24, z_offset=0.02) -> VGroup:
    V = np.asarray(V, dtype=np.float64)
    return VGroup(*[
        Text(str(i), font_size=font_size).move_to(p + np.array([0, 0, z_offset]))
        for i, p in enumerate(V)
    ])


def hinge_lines(V: NDArray[np.float64], edges, color=None) -> VGroup:
    V = np.asarray(V, dtype=np.float64)
    ls = []
    for (i, j) in edges:
        kwargs = {}
        if color is not None:
            kwargs["color"] = color
        ls.append(Line3D(V[int(i)], V[int(j)], **kwargs))
    return VGroup(*ls)
