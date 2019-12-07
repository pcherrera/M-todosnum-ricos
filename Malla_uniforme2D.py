import numpy as np
import GMSH_TRIANGULATION as msh

def shift(vector,k):
    return vector+k*np.array([1,1,1])

puntos = np.array([[0,0],[0,1],[1,1],[1,0]])
N = 4
x = np.linspace(0,1,N)
y = np.linspace(0,1,N)
malla=np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])
X = malla
# Elementos
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
print(X)
print(Ele)
print(len(Ele))
Th = msh.triangulation2DFEM(X, Ele)
# plot the mesh
msh.plotgmsh(Th.Coordinates, Th.Elements)
# Plot of the boundary
#msh.plotboundarygmsh(Th.Coordinates,Th.BoundaryEdges)
