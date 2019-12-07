import numpy as np
import scipy as sc
import generador as gen
class Malla1D:
    coordenadas = np.array([])
    numeracion_global = np.array([])
    elementos = np.array([])
    dirichlet = np.array([])
    neumann = np.array([])
    Ind = np.array([])
    N = 0
    def __init__(self,coordenadas,elementos,dirichlet,neumann):
        self.value = ''
        self.coordenadas = coordenadas
        self.elementos = elementos
        self.N = len(coordenadas)
        self.n_global = np.arange(0,self.N)
        self.dirichlet = dirichlet
        self.neumann = neumann
        self.Ind = np.setdiff1d(self.n_global,self.dirichlet)
    # Metodos
    def __set__(self, instance, value):
        self.value = value
    def __get__(self, instance, owner):
        return self.value

    def longitud(self,elemento): # Largo de un elemento
        a=elemento[0]
        b=elemento[1]
        p0 = self.coordenadas[a]
        p1 = self.coordenadas[b]
        return abs(p1-p0)






