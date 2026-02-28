# contains miscellanous functions and classes for testing syntax and python /numpy language knowledge 
import numpy as np
import manim 
from numpy.typing import NDArray

if __name__ == '__main__':
    a= np.array([1, 2, 3])
    a += 10
    print(a)

    A = np.ones((4,3))
    v = np.array([10, 20, 30])

    X = A + v
    print(X)