# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:55:06 2017

@author: bjohau
"""
import math
import numpy as np


def r_and_K(u, Lambda):
    r = np.array([
            [u[0,0] + 2.0 * u[0,0]**3 - u[1,0]**2 - 2.0 * Lambda],
            [3.0 * u[1,0] - 2.0 * u[0,0] * u[1,0] - Lambda]])
    
    K = np.array([
            [(1.0 + 6.0 * u[0,0]**2) ,    (-2.0 * u[1,0])],
            [      (-2.0 * u[1,0]) , (3.0 - 2.0 * u[0,0])]])
    
    return [r, K]






