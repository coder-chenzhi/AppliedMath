__author__ = 'chenzhi'

import matplotlib.pyplot as plt
import numpy as np

def draw_func(func, start=0, end=10):
    x = np.arange(start, end, 0.005)
    y = func(x)
    # print x
    # print y
    plt.plot(x, y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    # plt.show()

def main():
    draw_func()

if __name__ == "__main__":
    main()
