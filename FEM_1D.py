import numpy as np
import generador as gen
import matplotlib.pyplot as plt
import scipy.sparse.linalg as sci
class FEM1D:
    #ATRIBUTOS
    N = 1
    masa = np.array([])
    rigidez= np.array([])
    coordenadas = np.array([])
    elementos = np.array([])
    mreferencia = np.array([[1/3,1/6],[1/6,1/3]])
    kreferencia = np.array([[1,-1],[-1,1]])
    nodos_dirichlet = np.array([])
    nodos_independientes = np.array([])
    # Solucion
    solucion = np.array([])
    #SET
    def set_N(self,A): # |Dir|+|Ind|
        self.N = A
    def set_coordenadas(self,A):
        self.coordenadas = A
    def set_elementos(self,A):
        self.elementos = A
    def set_nodos_dirichlet(self,A):
        self.nodos_dirichlet = A
    def set_nodos_independientes(self,A):
        self.nodos_independientes = A
    #GET
    def get_coordenadas(self):
        return self.coordenadas
    def get_elementos(self,A):
        return self.elementos
    def get_nodos_dirichlet(self):
        return self.nodos_dirichlet
    def get_nodos_independientes(self):
        a = np.arange(0,self.N)
        return np.setdiff1d(a,self.nodos_dirichlet)

    # FUNCIONES
    def longitud(self,elemento): # Largo de un elemento
        a=elemento[0]
        b=elemento[1]
        p0 = self.coordenadas[a]
        p1 = self.coordenadas[b]
        return abs(p1-p0)
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
        #ALG GENERAL
        #for i in range(0,N):
            #for j in range(0,N):
                #A[i][j] = self.m(i,j)
        # MATRIZ TRIDIAGONA
        for i in range(0,N):
            A[i][i]= self.m(i,i)
        for i in range(0,N):
            A[i-1][i] = self.m(i-1,i)
            A[i][i-1] = self.m(i-1,i)
        #for i in range(0,N-2):
            #A[i+1][i] = self.m(i-1,i)

        return A
    def M_Ind(self): # Restriccion nodos independientes
        Free = self.nodos_independientes
        M = self.M()
        return M[Free,:][:,Free]
    def M_Dir(self): # Restriccion nodos dirichlet
        Ind = self.nodos_independientes
        Dir = self.nodos_dirichlet
        M = self.M()
        return M[Ind,:][:,Dir]


    # Matriz de Rigidez
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
        #for i in range(0,N):
            #for j in range(0,N):
                #A[i][j] = self.k(i,j)
        for i in range(0,N):
            A[i][i]= self.k(i,i)
        for i in range(0,N):
            A[i-1][i] = self.k(i-1,i)
            A[i][i-1] = self.k(i-1,i)

        return A
    def K_Ind(self): # Restriccion nodos Independientes
        Free = self.nodos_independientes
        K = self.K()
        return K[Free,:][:,Free]

    def K_Dir(self): # Restriccion Nodos Dirichlet
        Ind = self.nodos_independientes
        Dir = self.nodos_dirichlet
        K = self.K()
        return K[Ind,:][:,Dir]

    # Lado Derecho de la ecuación Ax = b

    def cos(self,x): # Funcion coseno
        return np.cos(x)

    def f_loc(self,K,i):
        if i in K :
            a = int(np.where(K == i)[0])
            p = self.coordenadas[i]
            longitud = self.longitud(K)
            return 0.5*longitud*2*self.cos(p)
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
        for i in range(0,N):
           A[i]= self.f(i)
        return A


    def b(self):
        #Generar vector (f,\phi_i)_i
        Ind = self.nodos_independientes
        Dir = self.nodos_dirichlet
        N = len(Ind)+len(Dir)
        A = np.zeros(N,dtype=np.float64)
        for i in range(0,N):
            A[i]= self.f(i)

        # Generar vector -u'(1)*\phi_i(1)
        a = np.zeros(N)
        a[N-1]=-np.sin(1)
        return A + a

    def b_ind(self):
        ind = self.nodos_independientes
        b = self.b()
        return b[ind]
    def g_dir(self):
        return 1

    # Cuadratura para (phi',cos')
    def g_loc(self,K,i):
        p0 = self.coordenadas[K[0]]
        p1 = self.coordenadas[K[1]]
        longitud = self.longitud(K)
        if i in K :
            a = int(np.where(K == i)[0])
            if a == 0:
                return (np.cos(p1)-np.cos(p0))/(p1-p0)
            else:
                return (np.cos(p0)-np.cos(p1))/(p1-p0)
        else:
            return 0

    def g(self,i):
        a = 0
        for k in self.elementos:
            a += self.g_loc(k,i)
        return a

    def G(self): # Vector F
        N = self.N
        A = np.zeros(N,dtype=np.float64)
        for i in range(0,N):
            A[i] = self.g(i)
        return A



    # Error L2 y H1
    def Error_L2(self):
        N = self.N
        U = self.solucion
        U_L2 = U@self.M()@U # norma L2 de U al cuadrado
        UF = np.sum(np.dot(self.F(),U)) # Producto escalar
        c = (1/4)*(2+np.sin(2)) #Valor exacto de la integral de cos^2
        return np.sqrt(U_L2-UF+c)

    def Error_H1(self):
        N = self.N
        U = self.solucion
        U_H1 = U@self.K()@U
        UF = 2*np.sum(np.dot(self.G(),U))
        c = (1/4)*(2-np.sin(2))
        return np.sqrt(U_H1-UF+c)
# PRUEBAS

FEM1D = FEM1D()
Error_L2 = np.array([])
Error_H1 = np.array([])

V = np.array([4,8,16,32,64,128,256])
Ratio = np.zeros(len(V))
for N in V:
    Ele = gen.uu_msh_ele1D(N)
    X = gen.uu_msh_points1D(N)
    Dir = np.array([0])
    # Iniciar FEM1D
    FEM1D.N = N
    FEM1D.coordenadas = X
    FEM1D.elementos = Ele
    FEM1D.nodos_dirichlet= Dir
    Ind = FEM1D.get_nodos_independientes()
    Dir = FEM1D.get_nodos_dirichlet()
    FEM1D.nodos_independientes= Ind

    # Relsolver sistema
    A = FEM1D.K_Ind()+FEM1D.M_Ind()
    b = FEM1D.b_ind()
    A_Ind = FEM1D.K_Ind()+FEM1D.M_Ind()
    A_dir = FEM1D.K_Dir()+FEM1D.M_Dir()
    F=b-np.transpose(A_dir)
    F = np.ravel(F)
    # Solucion
    U_Ind = np.linalg.solve(A_Ind,F)
    #U_Ind = sci.gmres(A_Ind,F)[0]
    U_Dir = np.array([1])
    m = len(Ind)+len(Dir)
    U = np.zeros(m)
    U[Ind] = U_Ind
    U[Dir] = U_Dir
    # Calcular error en L2
    FEM1D.solucion = U
    Error_L2 = np.append(Error_L2,FEM1D.Error_L2())
    Error_H1 = np.append(Error_H1,FEM1D.Error_H1())

for i in range(1,len(Error_L2)):
    r = Error_L2[i-1]/Error_L2[i]
    print('Error L2',r)
for i in range(1,len(Error_H1)):
    r = Error_H1[i-1]/Error_H1[i]
    print('Error H1',r)




