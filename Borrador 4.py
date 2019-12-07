import numpy as np
import GMSH_TRIANGULATION as gmsh
import generador as gen


class FEM2D(gmsh.triangulation2DFEM):
    # ATributos
    mreferencia = (1/24)*np.array([[2,1,1],[1,2,1],[1,1,2]])
    kxx = (1/2)*np.array([[1,-1,0],[-1,1,0],[0,0,0]])
    kyy = (1/2)*np.array([[1,0,-1],[0,0,0],[-1,0,1]])
    kxy = (1/2)*np.array([[1,0,-1],[-1,0,1],[0,0,0]])

    Masa = 0
    Rigidez = 0

    # MATRIZ DE MASA
    def M_K(self,i):
        N = len(self.Coordinates)
        B = np.zeros([N,N],dtype=np.float64) # matriz a rellenar

        K = self.Elements[i] # i-esimo triangulo
        Index = K
        # Auxiliares
        p0 = self.Coordinates[K[0],:]
        p1 = self.Coordinates[K[1],:]
        p2 = self.Coordinates[K[2],:]
        BK = np.array([p1-p0,p2-p0])
        det = abs(np.linalg.det(BK))
        # Rellenar matriz
        for a in range(0,3):
            for b in range(0,3):
                B[Index[a]][Index[b]] = det*self.mreferencia[a][b]
        return B

    def M(self):
        N = len(self.Coordinates)
        B = np.zeros([N,N],dtype=np.float64)
        r = len(self.Elements)
        for r in range(0,r):
            B += self.M_K(r)
        self.Masa = B
        return B

    # MATRIZ DE RIGIDEZ
    def K_K(self,i):
        N = len(self.Coordinates)
        B = np.zeros([N,N],dtype=np.float64) # matriz a rellenar

        K = self.Elements[i] # i-esimo triangulo
        Index = K
        # Auxiliares
        p0 = self.Coordinates[K[0],:]
        p1 = self.Coordinates[K[1],:]
        p2 = self.Coordinates[K[2],:]
        BK = np.array([p1-p0,p2-p0])
        BT = np.transpose(BK)
        C = np.linalg.inv(BT@BK)
        det = np.linalg.det(BK)
        s = self.kxy
        st = np.transpose(s)
        # Rellenar matriz
        for a in range(0,3):
            for b in range(0,3):
                B[Index[a]][Index[b]] = det*C[1][0]*(st[a][b]+s[a][b])\
                                        +det*C[0][0]*self.kxx[a][b] \
                                        +det*C[1][1]*self.kyy[a][b]
        return B

    def K(self):
        N = len(self.Coordinates)
        B = np.zeros([N,N],dtype=np.float64)
        r = len(self.Elements)
        for r in range(0,r):
            B += self.M_K(r)
        self.Masa = B
        return B
    # LADO DERECHO

    def f(self,x):
        return x[0]

    def g(self,x):
        return x[0]

    def G_L(self,i):
        N = len(self.Coordinates)
        L = self.Edges[i]
        B = np.zeros(N,dtype=np.float64)
        p0 = self.Coordinates[L[0]]
        p1 = self.Coordinates[L[1]]
        Largo = np.linalg.norm(p1-p0)
        for i in range(0,2):
            B[L[i]] = 0.5*Largo*self.g(0.5*(p0+p1)) # Regla de cuadratura
        return B

    def G(self):
        N = len(self.Coordinates)
        B = np.zeros(N,dtype=np.float64)
        r = np.where(self.BoundaryFlag == self.NeumannFlag)[0]
        Index = self.BoundaryEdgesIndex[r]

        for i in Index:
            B += self.G_L(i)
        return B


# FIN DE LA CLASE
#MALLA UNIFORME
#N = 32
#X = gen.uu_msh_points(N)
#Ele = gen.uu_msh_ele(N)

# OTRAS MALLAS
Th = gmsh.readmesh('Lshapeddomain5')
xy = [0,1] # solo queremos las coordenadas xy
# Iniciar FEM
FEM2D = FEM2D(Th.Coordinates[:][:,xy],Th.Elements,None,Th.BoundaryFlag,Th.DirichletFlag,Th.NeumannFlag)


