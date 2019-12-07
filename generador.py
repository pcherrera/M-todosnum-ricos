import numpy as np
import GMSH_TRIANGULATION as msh
def uu_msh_points(N): #unit uniform mesh
    x = np.linspace(0,1,N)
    y = np.linspace(0,1,N)
    malla=np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])
    return malla

def shift(vector,k):
    return vector+k*np.array([1,1,1])

def uu_msh_ele(N): #unit uniform elements
    puntos = np.array([[0,0],[0,1],[1,1],[1,0]])
    x = np.linspace(0,1,N)
    y = np.linspace(0,1,N)
    l = N**2 #nro de vertices
    a = np.array([0,1,N]) #Generadores
    b = np.array([1,N+1,N]) # Generadores
    Ele = np.array([a,b])
    for i in range(0,N-1):
        if i == 0 :
            for j in range(1,N-1):
                A = shift(a,j)
                Ele = np.vstack((Ele,A))
        else :
            for j in range(0,N-1):
                A = shift(a,N*i+j)
                Ele = np.vstack((Ele,A))

    for i in range(0,N-1):
        if i == 0 :
            for j in range(1,N-1):
                A = shift(b,j)
                Ele = np.vstack((Ele,A))
        else :
            for j in range(0,N-1):
                A = shift(b,N*i+j)
                Ele = np.vstack((Ele,A))
    return Ele


def uu_msh_points_dir(N):
    A = np.arange(0,N)
    return N*A

def uu_msh_edges_neu(N):
    A = np.array([0,1])
    for i in range(1,N-1):
        A = np.vstack((A,[i,i+1]))

    B = N*A
    C = N*(N-1)+A
    D = B+(N-1)
    A = np.vstack((A,C))
    A = np.vstack((A,D))
    return A

#def uu_msh_edges(Ele):



def uu_msh_ele1D(N):
    Ele = np.array([0,1])
    for i in range(1,N-1):
        Ele= np.vstack((Ele,[i,i+1]))
    return Ele
def uu_msh_points1D(N):
    return np.linspace(0,1,N)



#X = uu_msh_points(10)
#Ele = uu_msh_ele(10)
#msh.plotgmsh(X,Ele)
#msh.plotboundarygmsh(X,Ele)
#Edge =uu_msh_edges_neu(4)
#print(Edge)
#print(len(Edge))
#print(X[9])


