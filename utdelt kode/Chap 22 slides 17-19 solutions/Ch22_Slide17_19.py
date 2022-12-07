# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:55:06 2017

@author: bjohau
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import r_and_K as subs

#Parameter
a = 3

#Open output file
fp = open('Solution.txt','w')
#fp.print(fp,'Solution of 2x2 nonlinear problem\n\n')

#Guess initial solution value
u_vec = np.array([[0],[0]],float)
u_vec[0,0] = 0.8
u_vec[1,0] = 0.8  # Finds first solution
#u_vec[1,0] = 1.1   # Finds second solution
Lambda = 1.0

K_mat = np.zeros((2,2),float)
r_vec = np.zeros((2,1),float)

du_Vec = np.zeros((2,1))

nSteps = 10
for iLoop in range(nSteps):
    [r_vec, K_mat] = subs.r_and_K( u_vec, Lambda )
    
    du_vec = np.linalg.solve(K_mat,r_vec)
    
    u_vec = u_vec - du_vec
    
    print("i = {:2d}  x = [{:15.13f} {:15.13f}],  r = [{:15.13e} {:15.13e}]".format(iLoop, u_vec[0,0], u_vec[1,0], r_vec[0,0], r_vec[1,0]))
    





