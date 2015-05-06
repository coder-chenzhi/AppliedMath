__author__ = 'chenzhi'

"""
reference http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
"""

def calc_cost(x, y, cof):
    cost = 0.0
    for i in range(len(x)):
        tmp = 0.0
        for j in range(len(cof)):
            tmp += x[i] ** j * cof[j]
        cost += abs(tmp - y[i]) ** 2
    return cost

def calc_partial_derivatives(x, y, cof):
    partial_derivatives = []
    for k in  range(len(cof)):
        partial_derivative = 0.0
        for i in range(len(x)):
            tmp = 0.0
            for j in range(len(cof)):
                tmp += x[i] ** j * cof[j]
            partial_derivative += -2 * (y[i] - tmp) * (x[i] ** k)
        partial_derivatives.append(partial_derivative)
    return partial_derivatives

if __name__ == "__main__":
    print calc_cost([1,2,3],[1,2,3],[0, 1])
    print calc_partial_derivatives([1,2,3],[1,2,3],[0, 1, 1])
