## Elena Welch
## October 2022
## Implementation of fourth-order FDM scheme

import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt
from FDM_exactsolutions import u_exact1, u_exact2, du_exact1, du_exact2, Q_exact1, Q_exact2
from extrapolation_convergence import richardsons, convergence

L = 1
alphas = [0.2,0.4,0.75,3,5,7]
nodes = [2,4,8,16]
x_exact = np.linspace(0,L)
#################################################################
## FDM
## Case 1: U(0)=0, U(L)=0
def FDM_case1(n, L, alpha):
    deltax = L/n
    beta = alpha**2*deltax**2/12
    kappa = (2+10*beta)/(1-beta)
    K = np.zeros((n-1,n-1))
    f = np.zeros((n-1))
    f[-1] = -100
    K[-1,-1] = -kappa
    for i in range(n-2):
        K[i,i] = -kappa
        K[i,i+1] = 1
        K[i+1,i] = 1
    # print('K=',K)
    # print('f=',f)
    U = np.linalg.solve(K,f)
    U = np.concatenate(([0], U, [100]), axis=0)
    # print('U=', U)
    return U

## Case 2: U'(0)=0, U(L)=100
def FDM_case2(n, L, alpha):
    deltax = L/n
    beta = alpha**2*deltax**2/12
    kappa = (2+10*beta)/(1-beta)
    K = np.zeros((n,n))
    f = np.zeros((n))
    K[-1,-1]= -kappa
    f[-1] = -100
    for i in range(n-1):
        K[i,i] = -kappa
        K[0,0] = -kappa/2
        K[i,i+1] = 1
        K[i+1,i] = 1
    # print('K=',K)
    # print('f=',f)
    U = np.linalg.solve(K,f)
    U = np.concatenate((U, [100]), axis=0)
    # print('U=', U)
    return U

def grid_label(axis):
    for ax, col in zip(axis[0], nodes):
        ax.set_title("n = {}".format(col))

    for ax, row in zip(axis[:,0], alphas):
        a_label = r"$\alpha$ = {}".format(row)
        ax.set_ylabel(a_label, rotation=90, size='large')

## subplot set-ups
fig1,a1 = plt.subplots(len(alphas),len(nodes),squeeze=False,sharex=True,sharey=True,figsize=[8.5,11])
plt.suptitle("Case 1: FDM Temperature")
fig1.supxlabel("x")
fig1.supylabel("Temp (C)")
fig2,a2 = plt.subplots(len(alphas),len(nodes),squeeze=False,sharex=True,sharey=True,figsize=[8.5,11])
plt.suptitle("Case 2: FDM Temperature")
fig2.supxlabel("x")
fig2.supylabel("Temp (C)")

grid_label(a1)
grid_label(a2)

## loop through different values of alpha
## loop through different meshes for each alpha
u1_store = []  # store list for u1_exact at L/2
u2_store = []
U1_store = []
U2_store = []
for i in range(len(alphas)):
    print("alpha=",alphas[i])
    for j in range(len(nodes)):
        n = nodes[j]
        mesh = np.linspace(0,L,n+1)
        print("n =", n)
        U_FDM1 = FDM_case1(n, L, alphas[i])
        U_exact1 = u_exact1(mesh, alphas[i])
        error1 = np.abs((U_FDM1-U_exact1)/U_exact1)*100
        U_FDM2 = FDM_case2(n, L, alphas[i])
        U_exact2 = u_exact2(mesh, alphas[i])
        error2 = np.abs((U_FDM2-U_exact2)/U_exact2)*100
        if i == 0:
            U1_store.append(U_FDM1[n//2])
            U2_store.append(U_FDM2[n//2])

        print("Case 1")
        print(tabulate(zip(*[mesh, U_exact1, U_FDM1, error1]), headers=["x", "T_exact", "T_FDM", "Error"], tablefmt='latex'))
        print('\n')
        a1[i][j].plot(mesh,U_FDM1,label="FDM value",color='r',marker='x')
        a1[i][j].plot(x_exact,u_exact1(x_exact,alphas[i]), label='Exact',color='b')
        
        print("Case 2")
        print(tabulate(zip(*[mesh, U_exact2, U_FDM2, error2]), headers=["x", "T_exact", "T_FDM", "Error"], tablefmt='latex'))
        print('\n')
        a2[i][j].plot(mesh,U_FDM2,label="FDM value",color='g',marker='x')
        a2[i][j].plot(x_exact,u_exact2(x_exact,alphas[i]),label="Exact",color='k')

label_list = ['FDM value', 'Exact']

fig1.legend(labels=label_list,loc='upper left',borderaxespad=0.1)
fig1.savefig("FDM_fourthorder_case1.png")
fig2.legend(labels=label_list,loc='upper left',borderaxespad=0.1)
fig2.savefig("FDM_fourthorder_case2.png")
