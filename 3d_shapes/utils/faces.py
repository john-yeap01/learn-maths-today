import numpy as np
from numpy.typing import NDArray
from manim import Polygon


class Face:
    def __init__(self, name: str, V0: NDArray[np.float64], color, fill_opacity=0.25, stroke_width=2):
        self.name = str(name)
        self.V0 = np.asarray(V0, dtype=np.float64).copy()
        self.poly = Polygon(
            *[tuple(p) for p in self.V0],
            color=color,
            fill_opacity=fill_opacity,
            stroke_width=stroke_width,
        )
