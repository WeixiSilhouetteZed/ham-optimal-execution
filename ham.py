from sympy import symbols, poly, Sum, Indexed, sympify, Poly, Lambda, Derivative, Function, Subs, lambdify, DiracDelta, series
from sympy.functions import sign, Abs
import numpy as np
from scipy.special import gamma, factorial
import pandas as pd
import matplotlib.pyplot as plt

## initial Sympy setup


class HAM_Problem:
    def __init__(
        self, t_1: float, 
        t_2: float, 
        n: int, 
        delta: float,
        gamma: float,
        X: float):
        self.start_time = t_1
        self.end_time = t_2
        self.sub_num = n
        self.delta = delta 
        self.gamma = gamma
        self.X = X
        self.G_upper = np.fromfunction(self.G_ij_upper, (self.sub_num, self.sub_num))

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
        T = self.end_time - self.start_time
        N = self.sub_num
        return (1 / (1 - self.gamma) * (1 - self.gamma)) * (T / N) ** (2 - self.gamma) \
                * ((i - j + 1) ** (2 - self.gamma) - 2 * (i - j) ** (2 - self.gamma) +\
                    (i - j - 1) ** (2 - self.gamma)) 
    def G_ii_diagonal(self):
        T = self.end_time - self.start_time
        N = self.sub_num
        return (2 / ((1 - self.gamma) * (2 - self.gamma))) * (T / N) ** (2 - self.gamma)

    def create_ham_problem(self, h: float, M: int, init_type: str):
        self.M = M
        self.x, self.y, self.z, self.p, self.s, self.t = symbols("x y z p s t", real = True)
        self.v = symbols(f"v:{M}", real = True)
        self.v_sub = symbols(f"v_sub:{M}", real = True)
        self.d, self.g = symbols("delta gamma", real = True)
        f = Lambda(self.x, self.x ** self.d)
        f_dev = Lambda(self.x, self.d * self.x ** (self.d - 1))
        poly_p = Poly(self.v[::-1], self.p)
        poly_p_sub = Poly(self.v_sub[::-1], self.p)
        v_func = poly_p.expr
        v_sub_func = poly_p_sub.expr
        self.F_first = Lambda(self.p, f(v_func))
        self.F_second = Lambda(self.p, v_func * f_dev(v_sub_func))

        T, X, N = self.end_time - self.start_time, self.X, self.sub_num

        G_upper = np.fromfunction(self.G_ij_upper, (N, N))
        G_upper = np.nan_to_num(G_upper, nan = 0.0)
        G_diag = np.repeat(self.G_ii_diagonal(), N)
        G_diag = np.diag(G_diag)
        self.G = G_upper.T + G_upper + G_diag
        self.A = G_upper + G_diag / 2

        if init_type == "VWAP":
            v_guess = np.repeat(T * X / N, N) ## VWAP
        elif init_type == "GSS":
            v_guess = self.gss_guess(np.arange(0, 1 + 1 / (N + 1), 1 / (N + 1)))[1:-1] ## Linear

        self.v_guess = (v_guess / v_guess.sum()) * self.X
        update_array = np.zeros((N, M))
        update_array[:, 0] = self.v_guess
        self.result_array = self.ham_series_solve(update_array, h)
        self.solution = (self.result_array.sum(axis = 1) / self.result_array.sum()) * self.X
        self.expected_cost = self.expect_cost_compute(self.solution)

        # pd.Series(v_guess).plot(label = init_type)
        # pd.Series(self.solution).plot(label = "HAM")
        # plt.show()

        print(f"HAM solved for parameters delta = {self.delta}, gamma = {self.gamma}, guess = {init_type}")
        print(f"Theoretical cost = {self.expected_cost}")

    def gss_guess(self, t: np.array) -> np.array:
        T, X, N = self.end_time - self.start_time, self.X, self.sub_num
        gamma_1 = gamma((1 + self.gamma)/ 2)
        gamma_2 = gamma(1 + (self.gamma / 2))
        c = (X / N) / (np.sqrt(np.pi) * (T / 2) ** self.gamma * (gamma_1 / gamma_2))
        result = (c / (t * (T - t)) ** ((1 - self.gamma) / 2))
        return result

    def ham_series_solve(self, initial_array: np.array, h: float) -> np.array:
        update_dict = {self.p:0, self.d: self.delta}
        N = self.sub_num
        local_func = self.F_first(self.p)
        local_sub_func = self.F_second(self.p)
        for n in range(self.M - 1):
            for i in range(N):
                local_sum = 0
                for j in range(N):
                    local_g = self.G[i, j]
                    local_dict = update_dict.copy()
                    for idx in range(self.M):
                        local_dict[self.v[idx]] = initial_array[j, idx]
                        local_dict[self.v_sub[idx]] = initial_array[i, idx]
                    if j <= i:
                        local_f = local_func.subs(local_dict)
                    else:
                        local_f = local_sub_func.subs(local_dict)
                    local_sum += local_g * local_f
                initial_array[i, n + 1] = initial_array[i, n] + h * local_sum
            local_func = local_func.diff(self.p) / (n + 1)
            local_sub_func = local_sub_func.diff(self.p) / (n + 1)
        return initial_array

    def expect_cost_compute(self, v:np.array):
        return v @ self.A @ self.impact_f(v).T

def main():
    test_HAM = HAM_Problem(0, 1, 50, 0.5, 0.5, 0.1)
    test_HAM.create_ham_problem(-50.0, 6, "VWAP")
    print(test_HAM.solution)

if __name__ == "__main__":
    main()