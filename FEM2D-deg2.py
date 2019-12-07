import numpy as np
import generador as gen
import GMSH_TRIANGULATION as gmsh
import Cuadratura as Cuad
class FEM2D(gmsh.triangulation2DFEM):
    #Atributos
    masa = np.array([])
    rigidez= np.array([])
    ind_coord = np.array([])
    dir_coord = np.array([])
    neu_edges = np.array([])
    dir_val = np.array([])
    sol_val = np.array([])
    N = 0

    # OBJETOS DE REFERENCIA
    mreferencia = (1/360)*np.array([[6,-1,-1,-4,0,0],
                                   [-1,6,-1,0,-4,0],
                                   [-1,-1,6,0,0,-4],
                                   [-4,0,0,32,16,16],
                                    [0,-4,0,16,32,16],
                                    [0,0,-4,16,16,32]
                                    ])


    kxx = (1/6)*np.array([  [3,1,0,0,0,-4],
                            [1,3,0,0,0,-4],
                            [0,0,0,0,0,0],
                            [0,0,0,8,-8,0],
                            [0,0,0,-8,8,0],
                            [-4,-4,0,0,0,8]
                                    ])

    kyy = (1/6)*np.array([  [3,0,1,0,-4,0],
                            [0,0,0,0,0,0],
                            [1,0,3,0,-4,0],
                            [0,0,0,8,0,-8],
                            [-4,0,-4,0,8,0],
                            [0,0,0,-8,0,8]
                                    ])
    kxy = (1/6)*np.array([  [3,0,1,0,-4,0],
                            [1,0,-1,4,0,-4],
                            [0,0,0,0,0,0],
                            [0,0,4,4,-4,-4],
                            [0,0,-4,-4,4,4],
                            [-4,0,0,-4,4,4]
                                    ])

    # FUNCIONES NODALES
    triangle = np.array([[0,0],[0,1],[1,0]])

    # METODOS
    def is_dir(self,X):
        if X[0] == 0 and 0 <= X[1] <=1:
            return True
        else:
            return False

    def get_BoundaryFlags(self):
        B = np.array([])
        for i in self.BoundaryEdgesIndex:
            L = self.Edges[i]
            p0= self.Coordinates[L[0]]
            p1= self.Coordinates[L[1]]
            if (self.is_dir(p0) == True) and (self.is_dir(p1) == True):
                B = np.append(B,self.DirichletFlag)
            else:
                B = np.append(B,self.NeumannFlag)

        return B
    # MATRICES
    # Matriz de masa
    def M_K(self,i):
        N = len(self.Coordinates)
        L = len(self.Edges)
        B = np.zeros([N+L,N+L],dtype=np.float64) # matriz a rellenar

        K = self.Elements[i] # i-esimo triangulo
        Eleeges = self.EleEdges[i] + np.array([N,N,N])
        Index = np.concatenate((K,Eleeges))
        # Auxiliares
        p0 = self.Coordinates[K[0],:]
        p1 = self.Coordinates[K[1],:]
        p2 = self.Coordinates[K[2],:]
        BK = np.array([p1-p0,p2-p0])
        det = abs(np.linalg.det(BK))

        # Rellenar matriz
        for a in range(0,6):
            for b in range(0,6):
                B[Index[a]][Index[b]] = det*self.mreferencia[a][b]
        return B


    def M(self):
        N = len(self.Coordinates)
        L = len(self.Edges)
        K = len(self.Elements)
        B = np.zeros([N+L,N+L],dtype=np.float64)
        for i in range(0,K):
            B+= self.M_K(i)

        return B

    # Matriz de Rigidez

    def K_K(self,i):
        N = len(self.Coordinates)
        L = len(self.Edges)
        B = np.zeros([N+L,N+L],dtype=np.float64) # matriz a rellenar

        K = self.Elements[i] # i-esimo triangulo
        Eleeges = self.EleEdges[i] + np.array([N,N,N])
        Index = np.concatenate((K,Eleeges))

        #Auxiliares
        p0 = self.Coordinates[K[0],:]
        p1 = self.Coordinates[K[1],:]
        p2 = self.Coordinates[K[2],:]
        BK = np.array([p1-p0,p2-p0]).T
        BT = np.transpose(BK)
        C = np.linalg.inv(BT@BK)
        det = abs(np.linalg.det(BK))
        s = self.kxy
        st = np.transpose(s)


        # Rellenar matriz
        for a in range(0,6):
            for b in range(0,6):
                B[Index[a]][Index[b]] = det*C[1][0]*(st[a][b]+s[a][b])\
                                        +det*C[0][0]*self.kxx[a][b] \
                                        +det*C[1][1]*self.kyy[a][b]
        return B


    def K(self):
        N = len(self.Coordinates)
        L = len(self.Edges)
        K = len(self.Elements)
        B = np.zeros([N+L,N+L],dtype=np.float64)
        for i in range(0,K):
            B+= self.K_K(i)

        return B


    # Lado derecho
    def g_N(self,x):
        x1 = x[0]
        x2 = x[1]
        if (x1 == 1) and (0 < x2 <1):
            return -np.sin(1)
        else :
            return 0



    def F_K(self,i):
        N = len(self.Coordinates)
        L = len(self.Edges)

        K = self.Elements[i]
        Eleeges = self.EleEdges[i] + np.array([N,N,N])
        Index = np.concatenate((K,Eleeges))
        B = np.zeros(N+L,dtype=np.float64)

        #Auxiliares
        p0=self.Coordinates[K[0],:]
        p1=self.Coordinates[K[1],:]
        p2=self.Coordinates[K[2],:]
        Tri = np.array([p0,p1,p2])


        #detBk = abs(np.linalg.det(np.array([p2-p0,p1-p0]).T))

        for i in range(0,6):
            B[Index[i]] = Cuad.int_FK(Tri,i)
        return B

    def F(self):
        N = len(self.Coordinates)
        L = len(self.Edges)
        K = len(self.Elements)
        B = np.zeros(N+L,dtype=np.float64)
        for i in range(0,K):
            B+= self.F_K(i)
        return B

    def G_L(self,i):
        n = len(self.Coordinates)
        l = len(self.Edges)
        L = self.Edges[i]

        a = np.array([i])+n
        B = np.zeros(n+l,dtype=np.float64)

        Index = np.concatenate((L,a))

        p0 = self.Coordinates[L[0]]
        p1 = self.Coordinates[L[1]]

        largo = np.linalg.norm(p1-p0)
        for i in range(0,3):
            B[Index[i]] = 0.5*largo*self.g_N(0.5*(p0+p1))
        return B

    def G(self):
        n = len(self.Coordinates)
        l = len(self.Edges)
        B = np.zeros(n+l,dtype=np.float64)
        r = np.where(self.BoundaryFlag == self.NeumannFlag)[0]
        Index = self.BoundaryEdgesIndex[r]


        for i in Index:
            B += self.G_L(i)
        return B






#Fin de la clase

# Malla
N = 2
X = gen.uu_msh_points(N)
Ele = gen.uu_msh_ele(N)
# Iniciar FEM2D
FEM2D = FEM2D(X,Ele,None,None,1,2)
# Condicion de Neumann
FEM2D.BoundaryFlag = FEM2D.get_BoundaryFlags()


# Nodos y aristas de dirichlet e independientes
L = len(FEM2D.Edges)
dir_nodes = gen.uu_msh_points_dir(N)
r = np.where(FEM2D.BoundaryFlag == FEM2D.DirichletFlag)[0]
dir_edges = FEM2D.BoundaryEdgesIndex[r]+len(X)
Dir = np.concatenate((dir_nodes,dir_edges))
Ind = np.setdiff1d(np.arange(0,N**2+L),Dir)
# Valores dirichlet
U_dir = np.ones(len(Dir))


# Resolver sistema
F = FEM2D.F()
G = FEM2D.G()

b = FEM2D.F()+FEM2D.G()
A = FEM2D.M()+FEM2D.K()


U_ind = np.linalg.solve(A[Ind,:][:,Ind],b[Ind]-A[Ind,:][:,Dir]@U_dir)
U = np.zeros(len(Dir)+len(Ind))
U[Ind] = U_ind
U[Dir] = U_dir
#print(FEM2D.M())
# Establecer condiciones de borde
#print(FEM2D.F_K(1))
M = FEM2D.M()
print(M)
