import numpy as np
import matplotlib.pyplot as plt
Errores = np. array([
    [0.01468332039200049 , 0.08407772235039765],
    [0.002659054341209406 , 0.03533372573251674],
    [0.0005776312907008045 , 0.01642948571619776],
    [0.00013516898825348612 , 0.007943562996325],
    [3.2723925992344514e-05 , 0.0039080308678519325],
    [8.052413356727516e-06 , 0.0019385450449343892],
    [1.9973302635632965e-06 , 0.0009654611311165795]
    ])
print('Error L2')
for i in range(1,len(Errores)):
    print(np.log( Errores[i][0]/Errores[i-1][0] )/ np.log(2))
    #print(np.log(Errores[i][1])/np.log(Errores[i-1][1]))

print("Error H1")
for i in range(1,len(Errores)):
    print(np.log( Errores[i][1]/Errores[i-1][1] )/ np.log(2))
