import numpy as np
from Funciones1D import  Funciones1D
from scipy.sparse import spdiags

coordenadas = np.linspace(0,1,5)
FEM = Funciones1D()
FEM.set_coordenadas(coordenadas)
K = FEM.matriz_de_rigidez()
F = FEM.F()
C = FEM.matriz_principal()
U = np.linalg.solve(K+C,F)
FEM.set_solucion(U)
print("Matriz K \n",K)
print("Matriz M \n",C)
print("Vector F \n", F)
print("Solución U \n",U)
