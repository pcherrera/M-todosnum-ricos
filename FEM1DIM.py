import numpy as np
import generador as gen
import matplotlib.pyplot as plt
import scipy.sparse.linalg as sci
import quadpy as quad
class FEM1D:
    #ATRIBUTOS
    N = 1
    masa = np.array([])
    rigidez= np.array([])
    coordenadas = np.array([])
    elementos = np.array([])
    mreferencia = np.array([[1/3,1/6],[1/6,1/3]])
    kreferencia = np.array([[1,-1],[-1,1]])
    Dir = np.array([])
    Ind = np.array([])
    Neu = 0
    val_dir = 1
    val_neu = -np.sin(1)
    # Solucion
    solucion = np.array([])
    #SET

    # FUNCIONES
    def longitud(self,elemento): # Largo de un elemento
        a=elemento[0]
        b=elemento[1]
        p0 = self.coordenadas[a]
        p1 = self.coordenadas[b]
        return abs(p1-p0)
    def phi(self,L,i,u):
        x0 = L[0]
        x1 = L[1]
        t = (u-x0)/(x1-x0)
        if i == 0:
            return 1-t
        elif i == 1:
            return t

    def phi_u(self,L,i,u):
        x0 = L[0]
        x1 = L[1]
        t = (u-x0)/(x1-x0)
        if i == 0:
            return -1/(x1-x0)
        elif i == 1:
            return 1/(x1-x0)

    # Matriz de masa
    def loc_m(self,K,i,j):
        if i in K and j in K:
            a = int(np.where(K == i)[0])
            b = int(np.where(K == j)[0])
            longitud = self.longitud(K)
            return longitud*self.mreferencia[a][b]
        else:
            return 0

    def M(self):
        N = self.N
        Ele = self.elementos
        A = np.zeros([N,N],dtype=np.float64)
        for K in Ele:
            B = np.zeros([N,N],dtype=np.float64)
            for i in K:
                for j in K:
                    B[i][j]=self.loc_m(K,i,j)
            A +=B
        return A

    # Matriz de Rigidez
    def loc_k(self,K,i,j):
        if i in K and j in K:
            a = int(np.where(K == i)[0])
            b = int(np.where(K == j)[0])
            longitud = self.longitud(K)
            return (1/longitud)*self.kreferencia[a][b]
        else:
            return 0

    def K(self):
        N = self.N
        Ele = self.elementos
        A = np.zeros([N,N],dtype=np.float64)
        for K in Ele:
            B = np.zeros([N,N],dtype=np.float64)
            for i in K:
                for j in K:
                    B[i][j]=self.loc_k(K,i,j)
            A +=B
        return A

    def f_loc(self,K,i):
        if i in K :
            a = int(np.where(K == i)[0])
            p = self.coordenadas[i]
            longitud = self.longitud(K)
            return 0.5*longitud*2*np.cos(p)
        else:
            return 0
    def F(self):
        N = self.N
        Ele = self.elementos
        A = np.zeros(N,dtype=np.float64)
        for K in Ele:
            B = np.zeros(N,dtype=np.float64)
            for i in K:
                    B[i]=self.f_loc(K,i)
            A +=B
        C = np.zeros(N)
        C[self.Neu]=self.val_neu
        return A+C
    # Error
Error_L2 = np.array([])
Error_H1 = np.array([])
NN = np.array([4, 8, 16, 32,64,128,256,512])
FEM1D = FEM1D()
FEM1D.N = 4
for N in NN:

    X = np.linspace(0,1,N)
    Ele = gen.uu_msh_ele1D(N)
    # Iniciar Clase

    FEM1D.N = N
    FEM1D.coordenadas = X
    FEM1D.elementos = Ele
    A = FEM1D.M()+FEM1D.K()
    # Crear condiciones de borde
    Dir = np.array([0])
    Neu = np.array([N-1])
    Ind = np.setdiff1d(np.arange(0,N),Dir)
# Setear condiciones de borde
    FEM1D.Dir = Dir
    FEM1D.Ind = Ind
    FEM1D.Neu = Neu
    # Resolver sistema
    F = FEM1D.F()
    U_dir = np.array([1])
    U_ind = np.linalg.solve(A[Ind,:][:,Ind],F[Ind]-A[Ind,:][:,Dir]@U_dir)
    U = np.zeros(N)

    U[Ind] = U_ind
    U[Dir] = U_dir
    L2 = 0
    H1 = 0
    for j in range(0,len(FEM1D.elementos)):
        L = FEM1D.elementos[j]
        x0 = FEM1D.coordenadas[L[0]]
        x1 = FEM1D.coordenadas[L[1]]
        tri = np.array([x0,x1])
        funcion =  lambda x: (U[L[0]]*FEM1D.phi(tri,0,x)+U[L[1]]*FEM1D.phi(tri,1,x)-np.cos(x))**2
        funcion2 =  lambda x: (U[L[0]]*FEM1D.phi_u(tri,0,x)+U[L[1]]*FEM1D.phi_u(tri,1,x)+np.sin(x))**2
        Int_L = quad.line_segment.gauss_legendre(3).integrate(funcion,tri)
        Int_H = quad.line_segment.gauss_legendre(3).integrate(funcion2,tri)
        L2 += Int_L
        H1 += Int_H


    Error_L2 = np.append(Error_L2,L2)
    Error_H1 = np.append(Error_H1,H1)
Error_L2 = Error_L2**0.5
Error_H1 = Error_H1**0.5

r_L2 = np.zeros(len(Error_H1)-1)
r_H1 = np.zeros(len(Error_H1)-1)

for i in range(1,len(Error_L2)):
    r_L2[i-1] =  abs(np.log(Error_L2[i]/Error_L2[i-1])/np.log(2))
    r_H1[i-1] =  abs(np.log(Error_H1[i]/Error_H1[i-1])/np.log(2))


from tabulate import tabulate

hh = 1 / (NN - 1)
# Generamos la tabla de error
data = np.array([hh,
                 Error_L2 ,
                 np.concatenate((np.array(['-']), r_L2)),
                 Error_H1 ,
                 np.concatenate((np.array(['-']), r_H1))]
               ).T
head = ["h    ", "error L2 ", "r_L2", "error H1 ", "r_H1"]

print(tabulate(data, headers=head, tablefmt='fancy_grid', stralign='center'))
