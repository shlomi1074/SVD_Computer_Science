import numpy as np
from numpy import identity, diagonal, sort, unravel_index


def svd(A, tolerance=1.0e-8):
    At = A.T
    AtA = np.matmul(At, A)
    eigenValues, eigenVectors = calculate_eigens_jacobi(AtA, tolerance)

    index = eigenValues.argsort()[::-1]  # sort eigenValues in descending order and return array of matching index
    eigenValues = eigenValues[index]  # 'sort' eigenValues according to index array

    eigenVectors = eigenVectors[:, index]  # V matrix
    eigenValues = np.array([np.sqrt(i) if i >= 0 else 0 for i in eigenValues])

    eigenValuesInverse = np.divide(1, eigenValues, where=eigenValues != 0)

    eigenValuesMatrix = np.zeros([len(eigenValues), len(eigenValues)])  # S matrix
    np.fill_diagonal(eigenValuesMatrix, eigenValues)

    eigenValuesInverseMatrix = np.zeros([len(eigenValuesInverse), len(eigenValuesInverse)])
    np.fill_diagonal(eigenValuesInverseMatrix, eigenValuesInverse)
    '''A * eigenVectors * eigenValuesInverseMatrix'''
    U = np.matmul(np.matmul(A, eigenVectors), eigenValuesInverseMatrix)

    return U, eigenValuesMatrix, eigenVectors


def find_max_elem(A):
    upper_triangle = np.triu(A)  # set zeros in all values that not in upper triangle
    np.fill_diagonal(upper_triangle, 0.0)  # fill the diagonal with specified value - 0.0
    abs_upper_triangle = np.abs(upper_triangle)
    aMax = np.max(abs_upper_triangle)  # find the maximum value in numpy array

    # unravel_index get cell number and shape of input and return tuple (i,j) of element
    i_j = unravel_index(abs_upper_triangle.argmax(),
                        abs_upper_triangle.shape)  # argmax return the cell number of max element in matrix

    k = i_j.__getitem__(0)  # i = row
    l = i_j.__getitem__(1)  # j = column

    return aMax, k, l


def jacobi_eigens_rotate(AtA, maxRowIndex, maxColIndex, eigenVectors):
    maxRow = AtA[maxRowIndex][maxRowIndex]
    maxCol = AtA[maxColIndex][maxColIndex]
    maxRowCol = AtA[maxRowIndex][maxColIndex]

    aDiff = maxCol - maxRow
    if abs(maxRowCol) < abs(aDiff) * 1.0e-36:
        t = maxRowCol / aDiff
    else:
        phi = aDiff / (2.0 * maxRowCol)
        t = 1.0 / (abs(phi) + np.sqrt(pow(phi, 2) + 1.0))
        if phi < 0.0:
            t = -t

    c = 1.0 / np.sqrt(pow(t, 2) + 1.0)
    s = t * c
    tau = s / (1.0 + c)
    temp = maxRowCol

    '''updating M (for eigen values)'''
    AtA[maxRowIndex][maxColIndex] = 0
    AtA[maxRowIndex][maxRowIndex] = maxRow - t * temp
    AtA[maxColIndex][maxColIndex] = maxCol + t * temp

    for i in range(maxRowIndex):
        temp = AtA[i][maxRowIndex]
        elem = AtA[i][maxColIndex]
        tempValue = temp - s * (elem + tau * temp)
        AtA[i][maxRowIndex] = tempValue
        elem = AtA[i][maxColIndex]
        tempValue = elem + s * (temp - tau * elem)
        AtA[i][maxColIndex] = tempValue

    for i in range(maxRowIndex + 1, maxColIndex):
        temp = AtA[maxRowIndex][i]
        elem = AtA[i][maxColIndex]
        tempValue = temp - s * (elem + tau * temp)
        AtA[maxRowIndex][i] = tempValue
        elem = AtA[i][maxColIndex]
        tempValue = elem + s * (temp - tau * elem)
        AtA[i][maxColIndex] = tempValue

    for i in range(maxColIndex + 1, len(AtA)):
        temp = AtA[maxRowIndex][i]
        elem = AtA[maxColIndex][i]
        tempValue = temp - s * (elem + tau * temp)
        AtA[maxRowIndex][i] = tempValue
        elem = AtA[maxColIndex][i]
        tempValue = elem + s * (temp - tau * elem)
        AtA[maxColIndex][i] = tempValue

    '''updating eigenvectors'''
    for i in range(len(AtA)):
        tempRow = eigenVectors[i][maxRowIndex]
        tempCol = eigenVectors[i][maxColIndex]

        tempValue = tempRow - s * (tempCol + tau * tempRow)
        eigenVectors[i][maxRowIndex] = tempValue

        tempCol = eigenVectors[i][maxColIndex]
        tempValue = tempCol + s * (tempRow - tau * tempCol)
        eigenVectors[i][maxColIndex] = tempValue


def calculate_eigens_jacobi(AtA, tolerance=1.0e-8):
    nCols = len(AtA)
    eigenVectors = identity(nCols)
    numOfIterations = 2 * nCols
    for i in range(numOfIterations):
        aMax, maxValueRowIndex, maxValueColIndex = find_max_elem(AtA)
        if aMax < tolerance:
            break
        jacobi_eigens_rotate(AtA, maxValueRowIndex, maxValueColIndex, eigenVectors)

    return diagonal(AtA), eigenVectors
