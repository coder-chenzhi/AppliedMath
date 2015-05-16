__author__ = 'chenzhi'

"""
modified base on
http://scikit-learn.org/stable/auto_examples/mixture/plot_gmm.html
"""

import itertools

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn import mixture


def generate_data():
    n_samples = 500

    # Generate random sample, two components
    np.random.seed(0)

    # C is covariance matrix
    C = np.array([[0., -0.1], [1.7, .4]])

    # Translates slice objects to concatenation along the first axis
    # more details can be found in http://docs.scipy.org/doc/numpy/reference/generated/numpy.r_.html
    # you can do this by other more readable functions, like np.concatenate, more details can be found in
    # http://stackoverflow.com/questions/14468158/understanding-the-syntax-of-numpy-r-concatenation
    X = np.r_[np.dot(np.random.randn(n_samples, 2), C),
              .7 * np.random.randn(n_samples, 2) + np.array([-6, 3])]

    # np.dot(X, C) is dot product with covariance matrix C to do some shift
    X1 = np.r_[np.dot(np.random.randn(n_samples, 2), C)]

    # 0.7 is a multiplier to control variance, np.array([-6, 3]) is offset
    X2 = np.r_[.7 * np.random.randn(n_samples, 2) + np.array([-6, 3])]

    plt.scatter(X1[:, 0], X1[:, 1], .8, color="r")
    plt.scatter(X2[:, 0], X2[:, 1], .8, color="b")

    plt.show()

def main():
    # Number of samples per component
    n_samples = 500

    # Generate random sample, two components
    np.random.seed(0)
    C = np.array([[0., -0.1], [1.7, .4]])
    X = np.r_[np.dot(np.random.randn(n_samples, 2), C),
              .7 * np.random.randn(n_samples, 2) + np.array([-6, 3])]

    # Fit a mixture of Gaussians with EM using five components
    gmm = mixture.GMM(n_components=2, covariance_type='full')
    gmm.fit(X)
    print "GMM weights; ", gmm.weights_
    print "GMM means", gmm.means_
    print "GMM covariance", gmm.covars_

    # Fit a Dirichlet process mixture of Gaussians using five components
    dpgmm = mixture.DPGMM(n_components=5, covariance_type='full')
    dpgmm.fit(X)

    color_iter = itertools.cycle(['r', 'g', 'b', 'c', 'm'])

    for i, (clf, title) in enumerate([(gmm, 'GMM'),
                                      (dpgmm, 'Dirichlet Process GMM')]):
        splot = plt.subplot(2, 1, 1 + i)
        Y_ = clf.predict(X)
        for i, (mean, covar, color) in enumerate(zip(
                clf.means_, clf._get_covars(), color_iter)):
            v, w = linalg.eigh(covar)
            u = w[0] / linalg.norm(w[0])
            # as the DP will not use every component it has access to
            # unless it needs it, we shouldn't plot the redundant
            # components.
            if not np.any(Y_ == i):
                continue
            plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)

            # Plot an ellipse to show the Gaussian component
            angle = np.arctan(u[1] / u[0])
            angle = 180 * angle / np.pi  # convert to degrees
            ell = mpl.patches.Ellipse(mean, v[0], v[1], 180 + angle, color=color)
            ell.set_clip_box(splot.bbox)
            ell.set_alpha(0.5)
            splot.add_artist(ell)

        plt.xlim(-10, 10)
        plt.ylim(-3, 6)
        plt.xticks(())
        plt.yticks(())
        plt.title(title)

    plt.show()


if __name__ == "__main__":
    main()