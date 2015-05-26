__author__ = 'chenzhi'


def matcopy(A):
    return [a[:] for a in A]


def matinverse(AA, inplace=False):
    """
    Determines the inverse of a square matrix BB by Gauss-Jordan reduction.
    """
    n = len(AA)
    B = [[0] * n for i in range(n)]
    for i in range(n):
        B[i][i] = 1.0

    if not inplace:
        A = [row[:] for row in AA]
    else:
        A = AA
    try:
        for i in range(n):
            # Divide the ith row by A[i][i]
            m = 1.0 / A[i][i]
            for j in range(i, n):
                A[i][j] *= m  # # this is the same as dividing by A[i][i]
            for j in range(n):
                B[i][j] *= m

            #lower triangular elements.
            for k in range(i + 1, n):
                m = A[k][i]
                for j in range(i + 1, n):
                    A[k][j] -= m * A[i][j]
                for j in range(n):
                    B[k][j] -= m * B[i][j]

            #upper triangular elements.
            for k in range(0, i):
                m = A[k][i]
                for j in range(i + 1, n):
                    A[k][j] -= m * A[i][j]
                for j in range(n):
                    B[k][j] -= m * B[i][j]
        return B
    except:
        return None


def matprint(A, format="%8.3f"):
    # prints the matrix A using specified format
    m = len(A)
    try:
        n = len(A[0])
    except:
        n = 0  #
    if m:
        for i, row in enumerate(A):
            if n:
                for c in row:
                    print format % c,
                print
            else:
                print row
    print


def gausselim(AA, BB, pivoting=1, ztol=1.0e-8):
    """
     solves for X in AX = B.
     AA is a square matrix and B is a column vector(simple array).
     pivoting = 0.  no pivoting
              = 1.  partial row pivoting.
     Returns (code,R, A, X, B)
       codes: 0 - success,
              1 - near zero pivot encountered.Do not use the
                returned values!

    """
    A = [a[:] for a in AA]
    B = BB[:]

    size = len(A)
    X = [0.0] * size
    R = range(size)
    colperm = range(size)

    # Triangularization.
    for pivot in range(size):
        m = A[pivot][pivot]
        if pivoting == 1:
            absm = abs(m)
            for row in range(pivot + 1, size):
                testm = A[row][pivot]
                atestm = abs(testm)
                if atestm > absm:
                    m = (testm)
                    absm = abs(m)

                    # exchange pivot row with row row.
                    A[row], A[pivot] = A[pivot], A[row]
                    B[row], B[pivot] = B[pivot], B[row]

                    R[pivot], R[row] = R[row], R[pivot]

        if absm < ztol:
            return (1, R, A, X, B)  # missing row interchange vector in return.

        for row in range(pivot + 1, size):
            kmul = float(A[row][pivot]) / m
            # Apply rectangular rule.(Reduction)
            for col in range(size - 1, pivot, -1):
                A[row][col] = A[row][col] - kmul * A[pivot][col]
            B[row] = B[row] - kmul * B[pivot]
            A[row][pivot] = 0.0


            # Perform Back substitution.
    for row in range(size - 1, -1, -1):
        sum = B[row]
        for k in range(row + 1, size):
            sum -= (X[k] * A[row][k])
        # print row, sum,A[row][row]
        X[row] = sum / A[row][row]

    # print "Computed value of X = ",  X
    return (0, R, AA, X, B)


def SumSqr(X, Y, V, f, a):
    """
    Args
       X,Y - data vectors
       V   - data variance
       f   - function to minimize
       a   - fitting coefficients
    Auxiliary routine to compute the error sum of squares.
    """
    return sum([((y - f(x, a)) / v) ** 2 for x, y, v in zip(X, Y, V)])

    # previous lines.
    n = len(X)
    tot = 0.0
    for i in range(0, n):
        t = (Y[i] - f(X[i], a)) / V[i]
        tot += t * t
    return tot


def marquardt(X, Y, V, f, J, a, maxiter=20, minchange=1.0e-3, minlambdax=1.0e-6, debug=True, derh=0.04, dcode=1):
    """
    Args
      X, Y         Input:  data points (X_i, Y_i)
      V            Input:  data variance, must not have
                         zero values
      f            Input:  function with m unknown parameters
                         to be evaluated at point xi.
      J            Input:  d(f(X)) / d(a_k),
                         derivative of f() with respect to
                         the a_k parameter evaluated at
                         point x = xi, includes f() as
                         parameter for possible
                         numerical evaluation of derivative
      a            Input/Output: function parameters, will be updated by routine

      maxiter      maximum iterations
      minchange    Input:  minimum change between successive ChiSquares
      minlambda    Input:  minimum lambda value.
      debug        Input:  True if messages are to be printed.
      derh         Input:  numerical differentiation spacing.
      dcode        Input:  type of derivative approximation to apply.

    Return values
      flag
        0                 No error
        1                 Singular matrix
        2                 Maximum iterations exceeded
        3                 lambda becomes zero!
        4                 Singular covariance matrix !!
       besta              best fit parameters found by routine
       bestSS             best Sum of squares.
       Cov                covariance matrix
    """
    n = len(X)
    m = len(a)
    JtJ = [[0.0] * m for i in range(m)]

    bestSS = SS = SumSqr(X, Y, V, f, a)
    besta = a[:]
    Cov = None

    lambdax = 0.001
    iscomp = True
    ncount = 0
    flag = 0
    for p in range(1, maxiter + 1):
        if debug: print "marquardt(): iteration=", p
        if (iscomp):
            # Form JtJ matrix
            for k in range(m):
                for j in range(m):
                    tot = 0.0
                    for i in range(n):
                        tot += (J(X[i], f, a, k, derh, dcode) * J(X[i], f, a, j, derh, dcode)) / (V[i] * V[i])
                    JtJ[j][k] = JtJ[k][j] = tot

            if (lambdax == 0.0):
                if debug:
                    print "terminated, lambdax == 0!"
                break

            # Form RHS beta vector
            beta = [0] * m
            for k in range(m):
                beta[k] = sum([(y - f(x, a)) / (v * v) * J(x, f, a, k, derh, dcode) for x, y, v in zip(X, Y, V)])
        for j in range(m):
            JtJ[j][j] *= (1. + lambdax)


        # Solve for delta
        lastCov = matinverse(JtJ)
        if lastCov is not None:
            Cov = matcopy(lastCov)
        code, R, A, delta, b = gausselim(JtJ, beta)
        if code:
            flag = 4
            break
        totabsdelta = sum([abs(d) for d in delta])
        if debug:
            print "JtJ:"
            matprint(JtJ, "%5.3lf")
            print "beta = ", beta
            print "delta=", delta
            print "SS =", SS
            print "lambdax=", lambdax
            print "total abs delta=", totabsdelta

        if (code == 0):
            # Compute new parameters
            newa = [a[i] + delta[i] for i in range(m)]

            # and new sum of squares
            newSS = SumSqr(X, Y, V, f, newa)
            if debug: print "newSS = ", newSS
            # Update current parameter vector?
            if (newSS < bestSS):
                if debug: print "improved values found!"
                besta = newa[:]
                bestSS = newSS
                bestJtJ = matcopy(JtJ)
                a = newa[:]

                iscomp = True
                if debug:
                    print "new a:", a

                # Termination criteria
                if (SS - newSS < minchange):
                    ncount += 1
                    if (ncount == 2):
                        lambdax = 0.0
                        flag = 0
                        break
                else:
                    ncount = 0
                    lambdax = 0.4 * lambdax  # after Nash
                    if (lambdax < minlambdax):
                        flag = 3
                        break
                SS = newSS
            else:
                iscomp = False
                lambdax = 10.0 * lambdax
                ncount = 0
        else:
            flag = 1
            break

    if (flag == 0):
        if Cov is None: flag = 4
        if (p >= maxiter):
            flag = 2

    return flag, besta, bestSS, Cov


def J(xi, f, a, k, h=1.0e-5, code=2):
    """
    Args
       xi - evaluation point
       f  - function to minimize
       a  - fit parameters
       k  - index to active parameter.
       h  - small increment or spacing.
      code- approximation method
            0 - forward
            1 - backward
            2 - central
    Return
       numerical derivative of function f(xi, a) at point xi wrt to parameter a[k]
    """

    if code == 0:
        a[k] += h
        fxaph = f(xi, a)
        a[k] -= h  # restore ak
        fxa = f(xi, a)
        return (fxaph - fxa) / h

    elif code == 1:
        fxa = f(xi, a)
        a[k] -= h
        fxamh = f(xi, a)
        a[k] += h  # restore ak
        return (fxa - fxamh) / h

    elif code == 2:
        retval = 0.0

        a[k] += h
        retval = f(xi, a)
        a[k] -= 2 * h
        retval -= f(xi, a)
        a[k] += h  # restore ak
        return retval / (2 * h)


from math import exp


def f(xi, a):
    """
    Args
       xi  - evaluation point
       a[] - current fit parameters

    Desc
       Test function
    Return

       value of 2 parameter logistic function
    """
    return a[0] / (1. + a[1] * exp(xi * a[2]))


def f2(xi, a):
    return a[0] + a[1] / xi


if __name__ == "__main__":
    X = [1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.]
    Y = [5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948, 75.995, 91.972]
    V = [1.] * len(X)
    a = [0.7, 10., -0.4]  # initial estimates.

    # Sample Problem in
    X = range(1, 11)
    Y = [1.93, 7.13, 8.78, 9.69, 10.09, 10.42, 10.62, 10.71, 10.79, 11.13]
    a = [0, 0]
    derh = 1.0e-4
    dcode = 0  # use forward difference approximation.

    n = len(X)  # number of data points
    m = len(a)  # number of unknown function parameters
    maxiter = 50  # maximum iterations for matmarquardt()
    minchange = 0.1
    compcov = True

    print "\nInitial estimates:\n"
    for i in range(m):
        print "i = %d, %8.5lf\n" % (i, a[i])

    flag, a, minss, Cov = marquardt(X, Y, V, f2, J, a, maxiter, minchange, derh, dcode)

    print "Flag = %d, minSS = %f\n" % (flag, minss)

    for i in range(m):
        print "i = %d, %8.5lf\n" % (i, a[i])

    print "Last computed inverse"
    if flag == 0 and Cov != None:
        matprint(Cov, "%10.7lf ")
