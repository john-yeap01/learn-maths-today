# contains miscellanous functions and classes for testing syntax and python /numpy language knowledge 
import numpy as np
import manim 
from numpy.typing import NDArray


def stacker (v: list[float], w: list[float]):
    first  = np.array(v)
    print(first)
    print(f"list mode {v}")

    stacked = np.hstack((v,w))
    print(stacked)

    stacked_vertically = np.vstack((v,w))
    print( f"stacked vertically: {stacked_vertically}")


def my_matmul(v: list[float], w: list[float]):
    v = np.array(v)
    w = np.array(w)


    return v @ w


def broadcast(v: NDArray[np.float64]):
    return v + np.array([1,1,1])



if __name__ == '__main__':
    a: list[float] = [1,2,3]
    b: list[float] = [2,3,4]

    stacker(a,b)

    print(f"np array version of a: {np.array(a)}")


    c = np.array([[2,3,4], [3,4,5]])
    print(c.shape)

    d = c.ravel()
    print(f"d shape: {d.shape}")

    e = np.array([[3,4,5], [2,3,4]])
    print(e.shape)


    print(my_matmul(c,e.T))



    print(" broadcasting")
    print(broadcast(c))
    

    print(broadcast(a))
