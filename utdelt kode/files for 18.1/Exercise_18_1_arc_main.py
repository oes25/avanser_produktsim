# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:55:06 2017

@author: bjohau
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import K_and_q as subs


nSteps = 80
nVals = nSteps + 1
alpha = math.pi / 6
arc_length = 0.015
psi = 0
Lambda = 0  # Because lambda is reserved for small functions

psiVec = np.zeros([nVals, 1])
KVec = np.zeros([nVals, 1])
rVec = np.zeros([nVals, 1])
lambdaVec = np.zeros([nVals, 1])

lambdaVec[0] = Lambda
psiVec[0] = psi

d_psi_prev = 0.0
d_lambda_prev = 0.0

b_add_residual_solution = False

for iStep in range(0, nSteps):
    K, r = subs.K_func(psi, Lambda, alpha)
    # [K_num,r_num] = K_numerical(psi,Lambda,alpha)
    q = subs.q_func(psi, Lambda, alpha)
    v = q / K
    v_r = r / K

    # Length of augmented v vector [v 1]^T
    f = math.sqrt( v*v + 1)
    d_lambda = arc_length / f
    
    if ( (v * d_psi_prev + d_lambda_prev) < 0.0 ) :
        d_lambda *= -1.0
    
    d_psi = d_lambda * v
    
    Lambda = Lambda + d_lambda
    psi = psi + d_psi
    
#    if ( b_add_residual_solution ) :
#        psi = psi - v_r
    
    d_lambda_prev = d_lambda
    d_psi_prev = d_psi
    
    #Store for curve plotting
    lambdaVec[iStep+1] = Lambda
    psiVec[iStep+1] = psi
    KVec[iStep] = K
    rVec[iStep] = r

plt.plot(psiVec, lambdaVec)
#plt.show()
plt.plot(psiVec, rVec)
plt.grid()
plt.show()
#plt.legend('lambda')



