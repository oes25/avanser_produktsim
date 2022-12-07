# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 22:34:31 2018

@author: bjohau
"""

import math
import numpy as np
import matplotlib.pyplot as plt

a = 3

#file = open("Solution_from_python.txt","w")

nIter = 10

x = 1.0

for i in range(nIter):
    r = x**2 - a
    dr_dx = 2 * x
    dx = -r / dr_dx
    print("i = {:3d}  x = {:15.9e}   r = {:15.9e} dx = {:15.9e}".format(i, x, r, dx))
    dr_dx = 2 * x
    x = x + dx