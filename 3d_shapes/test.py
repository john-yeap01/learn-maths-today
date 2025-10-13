import numpy as np


def stacker (v: list[float], w: list[float]):
    first  = np.array(v)
    print(first)
    print(f"list mode {v}")



if __name__ == '__main__':
    a: list[float] = [1,2,3]
    b: list[float] = [2,3,4]

    stacker(a,b)