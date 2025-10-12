from manim import * 
import numpy
from typing import Tuple, Iterable

def unit(v: Iterable[float]) -> np.ndarray:
    x = np.asarray(v, dtype=float).reshape(3)

    # return np.array([1,2,3])
    return 2

class Test(ThreeDScene):
    def construct(self):
        test = unit([1,2,3])


