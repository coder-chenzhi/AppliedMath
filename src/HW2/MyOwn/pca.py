__author__ = 'chenzhi'

"""
reference:
http://www.cs.otago.ac.nz/cosc453/student_tutorials/principal_components.pdf
http://sebastianraschka.com/Articles/2014_pca_step_by_step.html

"""

import numpy as np
from matplotlib import pyplot as plt

def read_data(path):
    content = []
    data = []
    record = []
    with open(path) as f:
        content = f.readlines()
    for line in content:
        if line != "\n":
            record.extend([int(x) for x in line.strip("\n")])
        else:
            data.append(record)
            record = []
    return np.array(data).transpose()

if __name__ == "__main__":
    data = read_data("G:\Coding\Python\AppliedMath\src\HW2\MyOwn\\3")

    mean_vector = np.array([[np.mean(data[i, :])] for i in range(data.shape[0])])
    print('Mean Vector:\n', mean_vector)

    scatter_matrix = np.zeros((data.shape[0], data.shape[0]))
    for i in range(data.shape[1]):
        scatter_matrix += (data[:, i].reshape(data.shape[0], 1) -
                           mean_vector).dot((data[:, i].reshape(data.shape[0], 1) - mean_vector).T)
    print('Scatter Matrix:\n', scatter_matrix)

    # eigenvectors and eigenvalues for the from the scatter matrix
    eig_val_sc, eig_vec_sc = np.linalg.eig(scatter_matrix)
    # Make a list of (eigenvalue, eigenvector) tuples
    eig_pairs = [(np.abs(eig_val_sc[i]), eig_vec_sc[:,i]) for i in range(len(eig_val_sc))]

    # Sort the (eigenvalue, eigenvector) tuples from high to low
    # eig_pairs.sort()
    # eig_pairs.reverse()

    # Visually confirm that the list is correctly sorted by decreasing eigenvalues
    # for i in eig_pairs:
    #     print(i[0])
    matrix_w = np.hstack((eig_pairs[0][1].reshape(data.shape[0], 1), eig_pairs[1][1].reshape(data.shape[0], 1)))
    print matrix_w.shape
    print('Matrix W:\n', matrix_w)

    # Transforming the samples onto the new subspace
    transformed = matrix_w.T.dot(data)
    print transformed
    print transformed.shape
    plt.plot(transformed[0, :], transformed[1, :], 'o', markersize=3, color='green')
    plt.show()