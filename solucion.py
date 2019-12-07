import numpy as np
from Malla1D import Malla1D
import generador as gen
class solucion(Malla1D):
    sol = np.array([])
    val_dirichlet = np.array([])
    val_neumann = np.array([])
    K = np.array([]) # Rigidez
    M = np.array([]) # Masa
    mreferencia = np.array([[1/3,1/6],[1/6,1/3]])
    kreferencia = np.array([[1,-1],[-1,1]])


    # Matriz de masa
    def loc_m(self,K,i,j):
        if i in K and j in K:
            a = int(np.where(K == i)[0])
            b = int(np.where(K == j)[0])
            longitud = self.longitud(K)
            return longitud*self.mreferencia[a][b]
        else:
            return 0
    def m(self,i,j):
        A = 0
        for k in self.elementos:
            A += self.loc_m(k,i,j)
        return A
    def M(self):
        N = self.N
        A = np.zeros([N,N],dtype=np.float64)
        for i in range(0,N):
            for j in range(0,N):
                A[i][j] = self.m(i,j)
        return A

    # Matriz de Rigidez K
    def loc_k(self,K,i,j):
        if i in K and j in K:
            a = int(np.where(K == i)[0])
            b = int(np.where(K == j)[0])
            longitud = self.longitud(K)
            return (1/longitud)*self.kreferencia[a][b]
        else:
            return 0
    def k(self,i,j):
        A = 0
        for k in self.elementos:
            A += self.loc_k(k,i,j)
        return A
    def K(self):
        N = self.N
        A = np.zeros([N,N],dtype=np.float64)
        for i in range(0,N):
            for j in range(0,N):
                A[i][j] = self.k(i,j)
        return A
    #Lado derecho de Ax=b
    def f_loc(self,K,i): #(f,N_i)_K
        if i in K :
            p = self.coordenadas[i]
            a = int(np.where(K == i)[0])
            longitud = self.longitud(K)
            return 0.5*longitud*2*np.cos(p)
        else:
            return 0
    def f(self,i): # integral (f,\phi_i)
        a = 0
        for k in self.elementos:
            a += self.f_loc(k,i)
        return a

    def F(self): # Vector F
        N = self.N
        A = np.zeros(N,dtype=np.float64)
        B = np.zeros(N,dtype=np.float64)
        neu = self.neumann
        B[neu]=self.val_neumann
        for i in range(0,N):
           A[i]= self.f(i)
        return A+B


# Calcular solucion



