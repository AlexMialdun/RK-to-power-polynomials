# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 01:33:22 2018

@author: AlexMialdun
"""

import numpy as np
import numpy.polynomial.polynomial as poly
from scipy.optimize import minimize
from scipy.linalg import pascal


def polyvalRK(A, xdata):
    """
    Calculates a polynomial of Redlich-Kister type.
    Coefficients are defined by the vector 'A'.
    'xdata' vector is varying within [0 1].
    """
    xdata = np.array(xdata)
    Sum = np.zeros_like(xdata)
    for i in range(len(A)):
        Sum += A[i]*(1 - 2*xdata)**i
    return xdata*(1 - xdata)*Sum


def fitCoeffsRK(x, xdata, ydata):
    xdata = np.array(xdata)
    ydata = np.array(ydata)
    diff = ydata - polyvalRK(x, xdata)
    return np.sum(diff**2)


def polyfitRK(xdata, ydata, n):
    """
    Calculates a vector 'x' of coefficients of 
    a polynomial of Redlich-Kister type of (n-1)-th order.
    'xdata' vector is varying within [0 1].
    """
    x0 = np.ones(n)
    minRes = minimize(fitCoeffsRK, x0, args=(xdata, ydata), 
                      method='nelder-mead', 
                      options={'xtol': 1e-14, 'maxiter': 4000, 'disp': True})
    return minRes.x


def polyConvRK2power(A):
    """
    Convert full polynomial of Redlich-Kister type 
    into usual power series polynomial.
    """
    A = np.array(A)
    N = len(A)
    k = np.mgrid[0:N, 0:N][1]
    B = (-2)**k * pascal(N, kind='lower')
    Z = np.zeros((B.shape[0], 1))
    B = np.hstack((Z, B, Z)) - np.hstack((Z, Z, B))
    C = np.dot(A, B)
    return C


def polyConvRK2powerReduc(A):
    """
    Convert reduced polynomial of Redlich-Kister type 
    (without x*(1-x) prefix) into usual power series polynomial.
    """
    A = np.array(A)
    N = len(A)
    k = np.mgrid[0:N, 0:N][1]
    B = (-2)**k * pascal(N, kind='lower')
    C = np.dot(A, B)
    return C


def findExtremaRK(A):
    """
    Look for positions of extrema for the 
    Redlich-Kister polynomial given by 
    coefficients A. 
    """
    C = polyConvRK2power(A)
    polyRoots = poly.polyroots(poly.polyder(C))
    polyRoots = polyRoots[np.isreal(polyRoots)]
    #polyRoots = polyRoots[np.isclose(polyRoots.imag, 0)]
    Xi = polyRoots[(polyRoots.real>0) & (polyRoots.real<1)]
    return Xi.real
