
import math
import numpy as np

class OneDofSlantedTrussStructure:

    def __init__(self, LH, LV, K1, F0):
        self.LH = LH   # Horizontal length
        self.LV = LV   # Vertical height
        self.K1 = K1  # Stiffness of slanted truss
        self.F0 = F0  # Nominal load, to be scaled with lambda

        self.L0 = math.sqrt(LH**2 + LV**2) #Inital slanted length


    def get_vecQ(self,Lambda):
        vecQ = np.array([Lambda * self.F0])
        return vecQ

    def get_vecq(self,Lambda):
        vecq = np.array([self.F0])
        return vecq

    def get_vecInternalForces(self, u, Lambda):

        # Force from slanted part
        LV_d = self.LV + u[0]
        L_d = math.sqrt(self.LH**2 + LV_d**2)
        N1_slanted = (L_d - self.L0) * self.K1
        sinAlpha = LV_d / L_d
        N1_vertical = N1_slanted * sinAlpha

        vecInternalForce = np.array([-N1_vertical])
        return vecInternalForce

    def get_residual(self, u, Lambda):
        vecQ = self.get_vecQ(Lambda)
        vecInternal = self.get_vecInternalForces(u,Lambda)
        return (vecQ + vecInternal)

    def get_Kmat(self, u):
        Kmat = np.zeros((1,1))

        # Slanted element, material stiffness

        LV_d = self.LV + u[0]
        L_d = math.sqrt(self.LH**2 + LV_d**2)
        N1_slanted = (L_d - self.L0) * self.K1
        sinAlpha = LV_d / L_d
        cosAlpha = self.LH / L_d

        Kmat_vert = sinAlpha * self.K1 * sinAlpha
        Kmat[0,0] += Kmat_vert

        # Geometric stiffness
        N1_vertical = N1_slanted * sinAlpha

        Kgeo_vert = cosAlpha * N1_slanted * cosAlpha / L_d
        Kmat[0,0] += Kgeo_vert
        return Kmat




