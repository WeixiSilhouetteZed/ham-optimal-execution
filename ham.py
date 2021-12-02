from sympy import symbols, poly, Sum, Indexed, sympify, Poly, Lambda, Derivative, Function, Subs, lambdify, DiracDelta, series
from sympy.functions import sign, Abs
import numpy as np
from scipy.special import gamma, factorial
import pandas as pd

## initial Sympy setup
x, y, z, p, s, t = symbols("x y z p s t", real = True)


class HAM_Problem:
    def __init__(
        self, t_1: float, 
        t_2: float, 
        n: int, 
        v_0: np.array, 
        delta: float,
        gamma: float,
        X: float):
        self.start_time = t_1
        self.end_time = t_2
        self.sub_num = n
        self.v_0 = v_0 ## initial guess provided by bench mark strategy
        self.delta = delta 
        self.gamma = gamma
        self.X = X
        self.G_upper = np.fromfunction(self.G_ij_upper, (N, N))

    def time_split(self) -> np.array:
        t_1, t_2, n = self.start_time, self.end_time, self.sub_num
        N = (t_2 - t_1) / n
        return np.arange(t_1, t_2 + N, N)

    def F_ij(v: np.array, i: int, j: int) -> float:
        pass

    def impact_f(self, v: float) -> float: ## assume power-law impact
        return np.sign(v) * np.abs(v) ** self.delta
    
    def impact_f_dev(self, v_s, v_t):
        return v_s * self.delta * np.abs(v_t) ** (self.delta - 1)

    def G_ij_upper(self, i, j):
        T = self.t2 - self.t1
        N = self.sub_num
        return (1 / (1 - self.gamma) * (1 - self.gamma)) * (T / N) ** (2 - self.gamma) \
                * ((i - j + 1) ** (2 - self.gamma) - 2 * (i - j) ** (2 - self.gamma) +\
                    (i - j - 1) ** (2 - self.gamma)) 
    def G_ii_diagonal(self):
        T = self.t2 - self.t1
        N = self.sub_num
        return (2 / ((1 - self.gamma) * (2 - self.gamma))) * (T / N) ** (2 - self.gamma)

    def create_ham_problem(h: float, M: int):
        if 


def main():
    sample_ham = HAM_Problem(0, 100, 100, np.zeros(100), 0.5, 0.5)
    print(sample_ham.time_split())
    print(sample_ham.f(0.3))
    print(sample_ham.f_dev(0.4, 0.3))

if __name__ == "__main__":
    main()