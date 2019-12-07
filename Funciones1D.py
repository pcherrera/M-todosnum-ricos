from Malla1D import  Malla1D
import numpy as np
from scipy.sparse import spdiags
import quadpy as quad
import FEM1DIM as FEM
class fun:
    def phi(L,i,u):
        x0 = L[0]
        x1 = L[1]
        t = (u-x0)/(x1-x0)
        if i == 0:
            return 1-t
        if i == 1:
            return t


