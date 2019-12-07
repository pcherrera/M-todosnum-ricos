import numpy as np
import generador as gen
import quadpy as quad

class FEM2D:
    #Atributos
    masa = np.array([])
    rigidez= np.array([])
    coordenadas = np.array([])
    elementos = np.array([])
    ind_coord = np.array([])
    dir_coord = np.array([])
    neu_edges = np.array([])
    dir_val = np.array([])
    sol_val = np.array([])
    N = 0
    mreferencia = (1/24)*np.array([[2,1,1],[1,2,1],[1,1,2]])
    kxx = (1/2)*np.array([[1,-1,0],[-1,1,0],[0,0,0]])
    kyy = (1/2)*np.array([[1,0,-1],[0,0,0],[-1,0,1]])
    kxy = (1/2)*np.array([[1,0,-1],[-1,0,1],[0,0,0]])


    # Funciones
    def phi(self,Points,i,v):
        x0 = Points[0]
        x1 = Points[1]
        x2 = Points[2]
        det = (x1[0]-x0[0])*(x2[1]-x0[1])-(x1[1]-x0[1])*(x2[0]-x0[0])
        x = (1/det)*((x2[1]-x0[1])*(v[0]-x0[0])-(x2[0]-x0[0])*(v[1]-x0[1]))
        y = (1/det)*(-(x2[1]-x0[1])*(v[0]-x0[0])+(x1[0]-x0[0])*(v[1]-x0[1]))

        if i == 0:
            return 1-x-y
        elif i == 1:
            return x
        elif i == 2:
            return y

    def Gphi_x(self,Points,i,v):
        x0 = Points[0]
        x1 = Points[1]
        x2 = Points[2]


        det = (x1[0]-x0[0])*(x2[1]-x0[1])-(x1[1]-x0[1])*(x2[0]-x0[0])

        if i == 0 :
            x = (1/det)*((x2[1]-x1[1])-(x1[1]-x0[1]))
            #y = (1/det)*(-(x2[0]-x0[0])+(x1[0]-x0[0]))
            return  x
        elif i == 1:
            x = (1/det)*((x2[1]-x1[1]))
            #y = (1/det)*(-(x2[0]-x0[0]))
            return  x
        elif i == 2:
            x = (1/det)*(-(x1[1]-x0[1]))
            #y = (1/det)*((x1[0]-x0[0]))
            return x

    def Gphi_y(self,Points,i,v):
        x0 = Points[0]
        x1 = Points[1]
        x2 = Points[2]


        det = (x1[0]-x0[0])*(x2[1]-x0[1])-(x1[1]-x0[1])*(x2[0]-x0[0])

        if i == 0 :
            #x = (1/det)*((x2[1]-x1[1])-(x1[1]-x0[1]))
            y = (1/det)*(-(x2[0]-x0[0])+(x1[0]-x0[0]))
            return  y
        elif i == 1:
            #x = (1/det)*((x2[1]-x1[1]))
            y = (1/det)*(-(x2[0]-x0[0]))
            return  y
        elif i == 2:
            #x = (1/det)*(-(x1[1]-x0[1]))
            y = (1/det)*((x1[0]-x0[0]))
            return y




    # Metodos
    def m_loc(self,K,i,j):
        p0 = self.coordenadas[K[0],:]
        p1 = self.coordenadas[K[1],:]
        p2 = self.coordenadas[K[2],:]
        BK = np.array([p1-p0,p2-p0])
        a = int(np.where(K == i)[0])
        b = int(np.where(K == j)[0])
        det = np.linalg.det(BK)
        return det*self.mreferencia[a][b]

    def k_loc(self,K,i,j):
        p0 = self.coordenadas[K[0],:]
        p1 = self.coordenadas[K[1],:]
        p2 = self.coordenadas[K[2],:]
        BK = np.array([p1-p0,p2-p0]).T
        BT = np.transpose(BK)
        a = int(np.where(K == i)[0])
        b = int(np.where(K == j)[0])
        C = np.linalg.inv(BT@BK)
        det = np.linalg.det(BK)
        s = self.kxy
        st = np.transpose(s)
        return det*C[1][0]*(st[a][b]+s[a][b])\
               +det*C[0][0]*self.kxx[a][b] \
                +det*C[1][1]*self.kyy[a][b]

    def K(self):
        N = self.N
        A = np.zeros([N,N],dtype=np.float64)
        Ele = self.elementos
        for K in Ele:
            B = np.zeros([N,N], dtype=np.float64)
            for i in K :
                for j in K :
                    B[i][j] = self.k_loc(K,i,j)
            A += B
        return A

    def M(self):
        N = self.N
        A = np.zeros([N,N],dtype=np.float64)
        Ele = self.elementos
        for K in Ele:
            B = np.zeros([N,N], dtype=np.float64)
            for i in K :
                for j in K :
                    B[i][j] = self.m_loc(K,i,j)
            A += B
        return A

    # Lado derecho

    def cos(self,x): # Funcion coseno
        return np.cos(x)
    def g_N(self,x):
        x1 = x[0]
        x2 = x[1]
        if (x1 == 1) and (0 < x2 <1):
            return -np.sin(1)
        else :
            return 0

    def f_loc(self,K,i):
        p0=self.coordenadas[K[0],:]
        p1=self.coordenadas[K[1],:]
        p2=self.coordenadas[K[2],:]
        if i in K :
            detBk = abs(np.linalg.det(np.array([p2-p0,p1-p0]).T))
            #return (detBk/6)*(np.cos((p0[0]+p1[0]+p2[0])/3))
            return 2*(detBk/18)*(np.cos(p0[0])+np.cos(p1[0])+np.cos(p2[0]))
        else:
            return 0

    def F(self): # Vector F
        N = self.N
        A = np.zeros(N,dtype=np.float64)
        Ele = self.elementos
        for K in Ele:
            B = np.zeros(N,dtype=np.float64)
            for i in K:
                B[i]= self.f_loc(K,i)
            A += B
        return A

    def g_loc(self,L,i):
        p0=self.coordenadas[L[0]]
        p1=self.coordenadas[L[1]]
        largo = np.linalg.norm(p1-p0)
        return self.g_N(0.5*(p0+p1))*0.5*largo

    def h_loc(self,K,i):
        p0 = self.coordenadas[K[0],:]
        p1 = self.coordenadas[K[1],:]
        p2 = self.coordenadas[K[2],:]
        m1 = (p0+p1)/2
        m2 = (p1+p2)/2
        m3 = (p2+p0)/2
        grad1 = np.array([-np.sin(m1[0]),0])
        grad2 = np.array([-np.sin(m2[0]),0])
        grad3 = np.array([-np.sin(m3[0]),0])

        BK = np.array([p1-p0,p2-p0]).T
        BT = np.transpose(BK)
        a = int(np.where(K == i)[0])
        C = np.linalg.inv(BT@BK)
        det = abs(np.linalg.det(BK))

        if a == 0 :
            return -2*det*(C[0][0]+C[0][1])*(1-np.sin(1))
        elif a == 1 :
            return 2*det*C[0][0]*(1-np.sin(1))
        elif a == 2 :
            return 2*det*C[0][1]*(1-np.sin(1))

    def H(self): # Vector F
        N = self.N
        A = np.zeros(N,dtype=np.float64)
        Ele = self.elementos
        for K in Ele:
            B = np.zeros(N,dtype=np.float64)
            for i in K:
                B[i]= self.h_loc(K,i)
            A += B
        return A



    def G(self): # Vector F
        N = self.N
        A = np.zeros(N,dtype=np.float64)
        Edge = self.neu_edges
        for L in Edge:
            B = np.zeros(N,dtype=np.float64)
            for i in L:
                B[i]= self.g_loc(L,i)
            A += B
        return A

#Fin de la clase

Error_L2 = np.array([])
Error_H1 = np.array([])
NN = np.array([2,4, 8,16,32])
# Iniciar FEM2D

FEM2D = FEM2D()
FEM2D.N = 2
for N in NN :
    # Generar malla
    #Th = msh.readmesh('batman2')
    #Ele = Th.Elements.astype(int)
    #xy = [0,1]
    #X =Th.Coordinates[:][:,xy]
    FEM2D.N = 2
    X = gen.uu_msh_points(N)
    Ele = gen.uu_msh_ele(N)
    # Generar condiciones de borde
    Dir = gen.uu_msh_points_dir(N)
    neu_edges = gen.uu_msh_edges_neu(N)
    U_Dir = np.ones(len(Dir))
    # Setear malla
    FEM2D.coordenadas = X
    FEM2D.elementos = Ele
    FEM2D.N = len(X)
    # Setear condiciones de borde
    FEM2D.dir_coord = Dir
    FEM2D.dir_val = U_Dir
    FEM2D.neu_edges = neu_edges
    FEM2D.ind_coord = np.setdiff1d(np.arange(FEM2D.N),FEM2D.dir_coord)
    # Fabricar solucion
    Ind = FEM2D.ind_coord
    A = FEM2D.M()+FEM2D.K()
    F = FEM2D.F()+FEM2D.G()
    U_Ind = np.linalg.solve(A[Ind,:][:,Ind],F[Ind]-A[Ind,:][:,Dir]@U_Dir)
    U = np.zeros(N**2)
    U[Ind] = U_Ind
    U[Dir] = U_Dir

    L2 = 0
    H1 = 0
    # Calcular Error
    for j in range(0,len(Ele)):
        L = FEM2D.elementos[j]
        x0 = FEM2D.coordenadas[L[0]]
        x1 = FEM2D.coordenadas[L[1]]
        x2 = FEM2D.coordenadas[L[2]]
        p = (x1+x2+x0)/3
        tri = np.array([x0,x1,x2])
        funcion =  lambda x: (U[L[0]]*FEM2D.phi(tri,0,x)
                              +U[L[1]]*FEM2D.phi(tri,1,x)
                              +U[L[2]]*FEM2D.phi(tri,2,x)
                              -np.cos(x[0]))**2




        Int_L = quad.triangle.strang_fix_cowper_07().integrate(funcion,tri)
        L2 += Int_L

    H1 = U@FEM2D.K()@U-np.sum(FEM2D.H())+quad.line_segment.gauss_legendre(2).integrate(lambda x: np.sin(x)**2,[0,1])

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

