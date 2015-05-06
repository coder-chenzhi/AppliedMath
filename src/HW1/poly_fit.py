__author__ = 'chenzhi'

import numpy as np
import matplotlib.pyplot as plt
import random
import time


def raw_data(num, start=0, end=10):
    random.seed(time.time())
    x = sorted([random.uniform(start, end) for _ in range(num)])
    y = np.sin(x)
    return x, y


def gaussian_noise(num):
    noise = np.random.normal(0, 0.1, num)
    # print noise
    return noise


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
    f = np.polyfit(x, y, 2)
    draw_func(np.poly1d(f), label="fit")
    draw_func(np.sin, label="origin")
    plt.grid(True)
    plt.plot(x, y, marker='o', ls='')
    axes = plt.gca()
    axes.set_ylim([-2, 2])
    plt.show()

if __name__ == "__main__":
    main()