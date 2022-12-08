# example Beam_models
# ----------------------------------------------------------------
# PURPOSE
#  Starting point for a couple of beam models
#

import math
import numpy as np
import matplotlib.pyplot as plt
import CorotBeam
import matplotlib.animation as anm
from copy import deepcopy
# ----- Topology -------------------------------------------------


def solveArchLength(model, archLength=0.02, max_steps=50, max_iter=30):
    num_dofs = model.get_num_dofs()
    uVec = np.zeros(num_dofs)
    res_Vec = np.zeros(num_dofs)


    # v = displacement state
    # w = incremental velocity vector
    # w = dv/dLambda
    # q = incremental load vector
    # q = dp/dLambda
    # t = tangent of equilibrium path
    # f = length of tangent before normalization (kind of)

    # Getting stiffness matrix
    K_mat = model.get_K_sys(uVec)
    Lambda = 0.0
    q_Vec = model.get_incremental_load(Lambda)
    print("q_vec når lambda 0", q_Vec)

    print("Current q-vector: ", q_Vec)

    #
    d_q = np.linalg.solve(K_mat, q_Vec)
    print("d_q", d_q)

    uVec = d_q * Lambda

    # Unit tangent along equilibrium path
    w = np.array([0.5, 0.5])  # for example
    f = math.sqrt(1 + w.T @ w)

    d_q_prev = np.zeros(num_dofs)

    v_global = 0
    for iStep in range(max_steps):

        # TODO: Implement this

        delta_v = 0
        #for iIter in range(max_iter):
        for iIter in range(1):
            # TODO: Implement this

            # Finn et uttrykk for delta_v



            # regn ut residual
            res_Vec = model.get_residual(Lambda, uVec)

            K_mat = model.get_K_sys(uVec)
            # Displacement som resultat av residual
            d_r = np.linalg.solve(K_mat, res_Vec)

            # Displacement som resultat av q
            d_q = np.linalg.solve(K_mat, q_Vec)



            #uVec += delta_uVec

            #K =
            #d_r = np.linalg.solve(K, res_Vec)

            # Hvis residual liten nok -> delta_v OK
            if (res_Vec.dot(res_Vec) < 1.0e-15):
                break

        v_global += delta_v
        uVec = v_global

        model.append_solution(Lambda, uVec)
        print(" ")

def solveNonlinLoadControl(model, load_steps=0.01, max_steps=100, max_iter=30):
    num_dofs = model.get_num_dofs()
    uVec   = np.zeros(num_dofs)
    d_uVec = np.zeros(num_dofs)
    Lambda = 0


    for iStep in range(max_steps):

        Lambda = load_steps * iStep

        #TODO: Implement this

        for iIter in range(max_iter):

            # TODO: Implement this

            res_Vec = model.get_residual(Lambda, uVec)
            if (res_Vec.dot(res_Vec) < 1.0e-15):
                break

        model.append_solution(Lambda, uVec)
        print(" ")

        model.append_solution(Lambda, uVec)
        print("Non-Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))



def solveLinearSteps(model, load_steps=0.01, max_steps=100):
    num_dofs = model.get_num_dofs()
    uVec = np.zeros(num_dofs)

    for iStep in range(max_steps):

        Lambda = load_steps * iStep

        q_Vec = model.get_incremental_load(Lambda)
        print("q_Vec", q_Vec)
        K_mat = model.get_K_sys(uVec)

        d_q = np.linalg.solve(K_mat, q_Vec)

        uVec = d_q * Lambda

        model.append_solution(Lambda, uVec)
        print("Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))