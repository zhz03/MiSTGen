
# -*- coding: utf-8 -*-
"""
minimum snap single dimension geneator

@author: zhaoliang
"""
from qpsolvers import solve_qp
from numpy import array, dot
# from quadprog import solve_qp # this is not good enough 

if __name__ == '__main__':
    M = array([[1., 2., 0.], [-8., 3., 2.], [0., 1., 1.]])
    P = dot(M.T, M)  # quick way to build a symmetric matrix
    q = dot(array([3., 2., 3.]), M).reshape((3,))
    G = array([[1., 2., 1.], [2., 0., 1.], [-1., 2., -1.]])
    G = None
    h = array([3., 2., -2.]).reshape((3,))
    h = None
    A = array([1., 1., 1.])
    b = array([1.])
    
    x = solve_qp(P, q, None, None, A, b)
    # x = solve_qp(P, q, G, h)
    print(f"QP solution: x = {x}")

