# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:55:06 2017

@author: bjohau
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import K_and_q as subs


nSteps = 10
nVals = nSteps + 1
alpha = math.pi / 6
d_lambda = 0.1 / nSteps
psi = 0
Lambda = 0  # Because lambda is reserved for small functions

psiVec = np.zeros([nVals, 1])
KVec = np.zeros([nVals, 1])
rVec = np.zeros([nVals, 1])
lambdaVec = np.zeros([nVals, 1])

lambdaVec[0] = Lambda
psiVec[0] = psi

for iLoop in range(0, nSteps):
    K, r = subs.K_func(psi, Lambda, alpha)
    # [K_num,r_num] = K_numerical(psi,Lambda,alpha)
    q = subs.q_func(psi, Lambda, alpha)
    v = q / K
    Lambda = Lambda + d_lambda
    psi = psi + d_lambda * v

    lambdaVec[iLoop+1] = Lambda
    psiVec[iLoop+1] = psi
    KVec[iLoop] = K
    rVec[iLoop] = r

plt.plot(psiVec, lambdaVec)
#plt.show()
plt.plot(psiVec, rVec)
plt.grid()
plt.show()
#plt.legend('lambda')



