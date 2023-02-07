## Elena Welch
## September 2022
## AERO 430
## Richardson's Extrapolation

import numpy as np

def richardsons(Q_dx):  # where Q_dx = Qdx, Qdx/2, Qdx/4, ...
    Qextr = (Q_dx[1]**2-Q_dx[0]*Q_dx[2])/(2*Q_dx[1]-Q_dx[0]-Q_dx[2])
    return Qextr

def convergence(Q_dx, Qextr, dx):
    beta = np.log((Qextr-Q_dx[0])/(Qextr-Q_dx[1]))/np.log(2)
    return beta 

