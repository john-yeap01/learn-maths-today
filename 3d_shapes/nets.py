from manim import ThreeDScene
import numpy as np
from typing import Tuple, Iterable

def unit(v: Iterable[float]) :
    x = np.asarray(v, dtype=float).reshape(3)

    if x.size != 3:
        raise ValueError("Expected 3 components")
    
    n = float(np.linalg.norm(x))
    if n == 0.0:
        raise ValueError("zero length vector")
    
    return (x/n).astype(np.float64, copy=False)

def _skew(v: Iterable[float]):
    x = np.asarray(v, dtype=float)
    if x.size != 3:
        raise ValueError("Expected 3 components")
    x = x.reshape(3)
    # normalised: np.ndarray[Tuple[np.Any], np.dtype[np.float64]]= unit(x)
    x_, y_, z_ = map(float, x)
    skew = np.array([[0, -z_, y_], 
                                    [z_, 0, -x_],
                                    [-y_, x_, 0]], dtype=np.float64)
    return skew



def rodrigues(axis: Iterable[float], theta: float):
    k = unit(axis)
    K = _skew(k)

    s = float(np.sin(theta))
    c = float(np.cos(theta))

    I = np.eye(3,dtype=np.float64)

    R = I + s*K + (1-c)*(K@K)

    return R



class Test(ThreeDScene):
    def construct(self) -> None:
        K = _skew([1, 2, 3])
        assert np.allclose(K.T, -K)