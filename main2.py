import matplotlib.pyplot as plt
import numpy as np
from FDM_functions import FDM
from FDM_exactsolutions import *
from extrapolation_convergence import richardsons, convergence
from tabulate import tabulate

nodes = [2,4,8,16,32]
alphas = [0.2,0.4,0.75,3,5,7]
L = 1  # length of beam
k = 1  # heat flux (?) coefficient
A = 1  # area of beam

## plot exacts
x_exact = np.linspace(0,L)
## subplot set-up
def grid_label(axis):
    for ax, col in zip(axis[0], nodes):
        ax.set_title("n = {}".format(col))

    for ax, row in zip(axis[:,0], alphas):
        a_label = r"$\alpha$ = {}".format(row)
        ax.set_ylabel(a_label, rotation=90, size='large')

## loop through different values of alpha
## loop through different meshes for each alpha
## loop through different orders
for order in [2,4,6,8]:
    fig1,a1 = plt.subplots(len(alphas),len(nodes),squeeze=False,sharex=True,sharey=True,figsize=[8.5,11])
    plt.suptitle("Case 1: FDM Temperature")
    fig1.supxlabel("x")
    fig1.supylabel("Temp (C)")
    fig2,a2 = plt.subplots(len(alphas),len(nodes),squeeze=False,sharex=True,sharey=True,figsize=[8.5,11])
    plt.suptitle("Case 2: FDM Temperature")
    fig2.supxlabel("x")
    fig2.supylabel("Temp (C)")
    fig4,a4 = plt.subplots()
    fig3,a3 = plt.subplots()
    for i in range(len(alphas)): 
        U1_store = []
        U2_store = []
        deltax = []
        error1_store = []
        error2_store = []
        print("alpha=",alphas[i])
        for j in range(len(nodes)):
            n = nodes[j]
            deltax.append(L/n)
            mesh = np.linspace(0,L,n+1)
            print("n =", n)
            T_FDM1 = FDM(n, L, alphas[i],1,order)
            U_exact1 = u_exact1(mesh, alphas[i])
            error1 = np.abs((T_FDM1-U_exact1)/U_exact1)*100
            T_FDM2 = FDM(n, L, alphas[i],2,order)
            U_exact2 = u_exact2(mesh, alphas[i])
            error2 = np.abs((T_FDM2-U_exact2)/U_exact2)*100
            U1_store.append(T_FDM1[n//2])
            U2_store.append(T_FDM2[n//2])
            error1_store.append(error1[n//2])
            error2_store.append(error2[n//2])
            grid_label(a1)
            grid_label(a2) 
            print("Case 1")
            print(tabulate(zip(*[mesh, U_exact1, T_FDM1, error1]), headers=["x", "T_exact", "T_FDM", "Error"], tablefmt='github'))
            print('\n')
            a1[i][j].plot(mesh,T_FDM1,label="FDM value",color='r',marker='x')
            a1[i][j].plot(x_exact,u_exact1(x_exact,alphas[i]), label='Exact',color='b')
            
            print("Case 2")
            print(tabulate(zip(*[mesh, U_exact2, T_FDM2, error2]), headers=["x", "T_exact", "T_FDM", "Error"], tablefmt='github'))
            print('\n')
            a2[i][j].plot(mesh,T_FDM2,label="FDM value",color='g',marker='x')
            a2[i][j].plot(x_exact,u_exact2(x_exact,alphas[i]),label="Exact",color='k')
        a3.plot(-np.log(deltax), np.log(error1_store), label=r'$\alpha =$'+str(alphas[i]))
        a4.plot(-np.log(deltax), np.log(error2_store), label=r'$\alpha =$'+str(alphas[i]))
        U1_extr = richardsons(U1_store)
        U2_extr = richardsons(U2_store)
        beta1 = convergence(U1_store, U1_extr, mesh)
        beta2 = convergence(U2_store, U2_extr, mesh)
        print("beta1 = ", beta1)
        print("beta2 = ", beta2)

    label_list = ['FDM value', 'Exact']
    fig3.legend()
    a3.set(title="Rate of Convergence, Case 1",xlabel="-log(dx)",ylabel="log(error)")
    fig3.savefig("ROC_case1_"+str(order)+"order.png")
    fig4.legend()
    a4.set(title="Rate of Convergence, Case 2",xlabel="-log(dx)",ylabel="log(error)")
    fig4.savefig("ROC_case2_"+str(order)+"order.png")
    fig1.legend(labels=label_list,loc='upper left',borderaxespad=0.1)
    fig1.savefig("FDM_"+str(order)+"order_case1.png")
    fig2.legend(labels=label_list,loc='upper left',borderaxespad=0.1)
    fig2.savefig("FDM_"+str(order)+"order_case2.png")



