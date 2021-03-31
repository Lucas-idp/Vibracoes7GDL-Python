# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 18:11:20 2021

@author: Lucas
"""
from __future__ import print_function
import numpy as np
import sympy 
import math
from scipy.linalg import eigh
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.animation import FuncAnimation
from time import sleep

### Modeling System  

m_roda = 7 # wheel mass [kg]
m_carro = 230 #chassi mass [kg]
g = 9.81 #gravity [m/s^2]
L0 = 2.385 # total length[m]
L1 = (3/5)*1.7 # distance from front to CG [m]
L2 = 1.7-L1 #distance from rear to CG [m]
L3 = 1.4 #lateral wheet to wheel distance [m]

Ix_carro = 47 #rotational inercia of the chassi around X [kgm^2]
Iz_carro = 20 #rotational inercia of the chassi around Y [kgm^2]

k2a = 150*9.81/0.01 #spring constants 
k3a = 150*9.81/0.01
k4a = 150*9.81/0.01
k5a = 150*9.81/0.01
k2b = k2a/8
k3b = k3a/8
k4b = 1.5*k2b
k5b = 1.5*k3b

#### Declaring Symbolic Variables

y1_dot, y2_dot, y3_dot, y4_dot, y5_dot, gama_dot, teta_dot, y1, y2, y3, y4, y5, gama, teta = sympy.symbols('y1_dot, y2_dot, y3_dot, y4_dot, y5_dot, gama_dot, teta_dot, y1, y2, y3, y4, y5, gama, teta')      
var = [y2, y3, y4, y5, y1, gama, teta]
var_dot = [y2_dot, y3_dot, y4_dot, y5_dot, y1_dot, gama_dot, teta_dot]

#### Defining Kinetic and Potential Energy equations for the system

Ec = m_carro*(y1_dot**2)/2 + Ix_carro*(teta_dot**2)/2 + Iz_carro*(gama_dot**2)/2 + m_roda*((y2_dot**2)/2 + (y3_dot**2)/2 + (y4_dot**2)/2 + (y5_dot**2)/2) 
Ep = ((y1 + (gama*L1) + teta*(L3/2) - y3)**2)*k3b/2 + ((y1 - (gama*L2) + teta*(L3/2)- y2)**2)*k2b/2  + ((y1 + (gama*L2) - teta*(L3/2) - y4)**2)*k4b/2 + ((y1 - (gama*L2) - teta*(L3/2) - y5)**2)*k5b/2  + (y2**2)*k2a/2 + (y3**2)*k3a/2 + (y4**2)*k4a/2 + (y5**2)*k5a/2 + m_roda*g*(y1 + y2 + y3 + y4 + y5)


##### Defining Lagragian as function of the two energies
L = Ec - Ep
L = sympy.simplify(L)

#### Derivatives with respect to displacements
dL_dy = []
#### Derivatives with respect to velocities
dL_dy_dot = []
#### Writing equations of motion for every degree of freedom

dL_dy_dot = []
graus = []

for i in range(len(var)):
    dL_dy.append(0)
    dL_dy_dot.append(0)
    graus.append(0)
    dL_dy[i] = sympy.diff(L, var[i])# Derivatives with respect to displacements
    dL_dy_dot[i] = sympy.diff(L, var_dot[i])# #### Derivatives with respect to velocities
    graus[i] = dL_dy_dot[i] + dL_dy[i] # Writing equations of motion for every degree of freedom
    print(graus[i],'\n')
    
### Assembly of rigidity matrix
K = []

for i in range(len(var)):
    k = []
    for j in range(len(var)):
        k.append(0)
    K.append(k)
    
print(K)

for i in range(len(var)):
    for j in range(len(var)):
       
       K[i][j] = graus[i].coeff(var[j],1)

for i in range(len(var)):
    for j in range(len(var)):
        if K[i][j] == K[j][i]:
           pass
        else:
           print('Error: Matrix is not symmetric')
           break
print('The Matrix is symmetric') 

### Assembly of mass matrix

M = []

for i in range(len(var)):
    m = []
    for j in range(len(var)):
        m.append(0)
    M.append(m)
    
print(M)

for i in range(len(var)):
    for j in range(len(var)):
       
       M[i][j] = graus[i].coeff(var_dot[j],1)
print(M)
for i in range(len(var)):
    for j in range(len(var)):
        if M[i][j] == M[j][i]:
           pass
        else:
           print('Error: Matrix is not symmetric')
           break
print('The Matrix is symmetric')


####3 Modal Analysis

##### Eigenvalues and Eigenvectors
for i in range(len(var)):
    for j in range(len(var)):
        K[i][j] = float(K[i][j])
        M[i][j] = float(M[i][j])
        
M = np.array(M)
K = np.array(K)        
[wn, B] = eigh(K,M)
wn = np.sqrt(wn**2)
print(wn)
print(B) 

freq_n = np.sqrt(((wn)**2)**(1/2))
print(freq_n)

### Assembly of modified mass matrix or modal mass matrix 
meff = []
for i in range(len(var)):
    meff.append(0)
    meff[i] = np.dot(np.dot(B[:,i].T,M),B[:,i])
   
meff = np.diag(meff)
for i in range(len(var)):
    for j in range(len(var)):
        if meff[i][j] == meff[j][i]:
           pass
        else:
           print('Error: Matrix is not symmetric')
           break
print('The Matrix is symmetric')

 #### initial conditions
x0 = [0,0.3,0,0,0,0,0]
v0 = [0,0,0,0,0,0,0]
y0 = np.dot(np.dot(B.T,M),x0)
y0dot = np.dot(np.dot(B.T,M),v0)
# total analysis time and time increments
t_inc0 = 0.005 #seconds
t_end0 = 10 #seconds   
y = []
t = np.arange(0, t_end0, t_inc0)
for r in range(len(var)):
    y.append(0)
    y[r] = y0[r]*np.cos(freq_n[r]*t) + y0dot[r]/(freq_n[r]*np.sin(freq_n[r]*t))
    
x = np.dot(B,y)
x1 = x[0][:]
x2 = x[1][:]
x3 = x[2][:]
x4 = x[3][:]
x5 = x[4][:]
x6 = x[5][:]
x7 = x[6][:]
for i in range(len(var)):
    plt.plot(t,x[i][:])
np.savetxt('7GDL_RespostaAnimBlue', np.transpose([x1,x2,x3,x4,x5,x6,x7]), fmt='%1.5f')





