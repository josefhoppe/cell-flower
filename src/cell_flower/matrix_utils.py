"""
Some helpers for linalg / matrix operations
"""

import numpy as np
import numpy.linalg as linalg


def normalize(matrix: np.ndarray) -> np.ndarray:
    interm = np.diag(np.diag(np.ones(matrix.shape).T @ matrix))
    return matrix @ linalg.inv(interm)

def neg_part(matrix: np.ndarray) -> np.ndarray:
    return -(np.minimum(matrix, np.zeros(matrix.shape)))

def pos_part(matrix: np.ndarray) -> np.ndarray:
    return np.maximum(matrix, np.zeros(matrix.shape))
