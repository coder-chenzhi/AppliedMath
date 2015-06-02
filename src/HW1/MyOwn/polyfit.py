__author__ = 'chenzhi'

"""
reference http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
"""

import numpy as np
import matplotlib.pyplot as plt
import random

def calc_poly_func(point, cof):
    v = 0
    for i in range(len(cof)):
        v += point ** i * cof[i]
    return v

def calc_cost(x, y, cof):
    cost = 0.0
    for i in range(len(x)):
        y_tmp = calc_poly_func(x[i], cof)
    try:
        cost += abs(y_tmp - y[i]) ** 2
    except Warning:
        pass
        # print "calc cost overflow", tmp, x, y, cof, cost
    return cost


def calc_partial_derivatives(x, y, cof):
    partial_derivatives = []
    for i in range(len(cof)):
        partial_derivative = 0.0
        for j in range(len(x)):
            tmp = calc_poly_func(x[j], cof)
            try:
                partial_derivative -= 2 * (y[j] - tmp) * (x[j] ** i)
                # print "partial_derivative", partial_derivative
            except Warning:
                pass
                # print "calc partial_derivative overflow x[i], y[i], k, tmp", x[j], y[j], i, tmp
        partial_derivatives.append(partial_derivative)
    return partial_derivatives


def poly_fit(x, y, deg, alpha=0.000001, rounds=1000, prec=0.01):
    initial_point_num = 1
    minimum_cost = float("inf")
    best_cof = []
    for _ in range(initial_point_num):
        cof = (np.random.random(deg) - 0.5).tolist()
        cur_round = 0
        while calc_cost(x, y, cof) > prec and cur_round < rounds:
            print cur_round, cof
            partial_derivatives = calc_partial_derivatives(x, y, cof)
            # print "partial_derivatives", partial_derivatives
            for i in range(len(cof)):
                cof[i] += alpha * partial_derivatives[i]
            cur_round += 1

        if calc_cost(x, y, cof) < minimum_cost:
            best_cof = cof
    return best_cof


def calc_func(cof, x):
    y = []
    for i in range(len(x)):
        value = 0.0
        for j in range(len(cof)):
            value += x[i] ** j * cof[j]
        y.append(value)
    return y


def raw_data(num, start=0, end=10):
    random.seed(17)
    x = sorted([random.uniform(start, end) for _ in range(num)])
    y = np.sin(x)
    return x, y


def gaussian_noise(num):
    noise = np.random.normal(0, 0.1, num)
    # print noise
    return noise

def draw_func_cof(cof, start=0, end=10, label=None):
    x = np.arange(start, end, 0.005)
    y = calc_func(cof, x)
    print "x", x
    print "y", y

    # max_y = max(y)
    # min_y = min(y)
    # ylim = max(abs(max_y), abs(min_y)) + 1

    if label is not None:
        plt.plot(x, y, label=label, )
    else:
        plt.plot(x, y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    # plt.show()


def draw_func(func, start=0, end=10, label=None):
    x = np.arange(start, end, 0.005)
    y = func(x)
    # print x
    # print y

    # max_y = max(y)
    # min_y = min(y)
    # ylim = max(abs(max_y), abs(min_y)) + 1

    if label is not None:
        plt.plot(x, y, label=label, )
    else:
        plt.plot(x, y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    # plt.show()


def main():
    num = 10
    x, y = raw_data(num)
    y += gaussian_noise(num)
    x = [1, 2, 3]
    y = [1, 2, 3]
    cof = poly_fit(x, y, 2, alpha=0.0001, rounds=10000, prec=0.01)
    print "cof", cof
    draw_func_cof(cof, label="fit")
    # draw_func(np.sin, label="origin")
    plt.grid(True)
    plt.plot(x, y, marker='o', ls='')
    axes = plt.gca()
    axes.set_ylim([-2, 2])
    plt.show()

if __name__ == "__main__":
    main()
