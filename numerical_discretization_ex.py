# projeto IC Levitação Acústica

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#teste 


#Dados

g = 9.81 # m/s²
sin_alpha = np.sin(np.pi/6) # rad
mu = 1.7e-1 # kg/m*s  viscosidade cinemática óleo SAE 10W30
b = 0.5 # m²
h = 2e-2 # m
m = 7 # kg

#Dados para a discretização da EDO 1° Ordem ñ homogenea

t = 0 # s
tf = 35# s
# o incremento de tempo é dt e deve ser uma vari-
# avel de entrada

dt = 0.005 # s

Vf = 0
V = 0 # m/s
mtx = np.zeros(((int((tf-t)/dt)),2))

for i in range(int((tf-t)/dt)): # no range do numero de passos

    mtx[i][0] = t
    mtx[i][1] = V

    Vf = V + (dt * ((g * sin_alpha) - (V * ((mu * b) / (h*m)))))
    t = t + dt
    V = Vf




plt.figure()
plt.plot(mtx[:,0],mtx[:,1])
plt.show()
