## Elena Welch
## September 2022
## AERO 430
## HW 2 

from FDM_fourthorder import U1_store, U2_store, nodes, L, error1, error2
from extrapolation_convergence import richardsons, convergence
from ROC import rateofconv
import matplotlib.pyplot as plt


dx = [L/n for n in nodes]
extr_value1 = richardsons(U1_store)
extr_value2 = richardsons(U2_store)
# beta1 = convergence(U1_store,extr_value1,dx)
# beta2 = convergence(U2_store,extr_value2,dx)
beta1 = rateofconv(error1, dx)
beta2 = rateofconv(error2,dx)
print("extrapolated temp, case 1 = ", extr_value1)
print("extrapolated temp, case 2 = ", extr_value2)
print("beta: case 1 = ", beta1)
print("beta: case 2 = ", beta2)

