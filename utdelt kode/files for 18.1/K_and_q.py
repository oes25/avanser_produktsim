# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:55:06 2017

@author: bjohau
"""
import math


def K_func(psi, Lambda, alpha):
    # K_val = partial(r)/partial(psi)
    # r_val = residual r

    r1 = (1/math.cos(alpha - psi))
    dr1 = - (1/math.cos(alpha - psi)) * math.tan(alpha - psi)

    a1 = (1/math.cos(alpha)) * (2 + Lambda * math.sin(psi)) - 2 / math.cos(alpha - psi)
    da1 = (1/math.cos(alpha)) * Lambda * math.cos(psi) + 2 * math.tan(alpha - psi)/ math.cos(alpha - psi)

    a2 = math.tan(alpha - psi)
    da2 = - (1/math.cos(alpha - psi)) ** 2

    a3 = Lambda * math.cos(psi) / math.cos(alpha)
    da3 = - Lambda * math.sin(psi) / math.cos(alpha)

    r2 = a1 * a2 - a3
    dr2 = da1 * a2 + a1 * da2 - da3

    r = r1 * r2
    dr = dr1 * r2 + r1 * dr2

    return [dr, r]


def q_func(psi, Lambda, alpha):
    # q = - partial(r)/partial(lambda)

    # sec = 1 / cos
    r1 = 1 / math.cos(alpha - psi)
    dr1 = 0

    a1 = (1 / math.cos(alpha)) * (2 + Lambda * math.sin(psi)) - 2 / math.cos(alpha - psi)
    da1 =  math.sin(psi)/ math.cos(alpha)

    a2 = math.tan(alpha - psi)
    da2 = 0

    a3 = Lambda * math.cos(psi) / math.cos(alpha)
    da3 = math.cos(psi) / math.cos(alpha)

    r2 = a1 * a2 - a3
    dr2 = da1 * a2 + a1 * da2 - da3

    dr = dr1 * r2 + r1 * dr2
    q_val = -dr

    return q_val




