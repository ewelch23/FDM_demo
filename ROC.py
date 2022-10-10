import numpy as np
import matplotlib.pyplot as plt

L = 1
u_exact1 = lambda x,alpha: 100/np.sinh(alpha * L) * np.sinh(alpha * x)
def FDM_case1(n, L, alpha):
    deltax = L/n
    kappa = 2+alpha**2*deltax**2
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

meshes = [2,4,8,16,32,64,128]
deltax = []
perror1 = []
U1_FDM = []
uexact = []
U1_ROC = []
u1_ROC = []
for i in range(len(meshes)):
    n = meshes[i]
    print("n=", n)
    deltax.append(L/n)
    U = FDM_case1(n,L,3)
    print("U = ", U)
    for j in range(1,n):
        u = u_exact1(j*L/n, 3)
    print("u=",u)
    print("##########")
    U1_FDM.append(U)
    uexact.append(u)
    perror1.append(np.abs((U-u)/u) * 100)

def rateofconv(errors,deltax):  # richardsons extrapolation: compute betas
    betas = []
    for i in range(1,len(errors)):
        beta = np.log(errors[i]/errors[i-1])/np.log(deltax[i]/deltax[i-1])
        betas.append(beta)
    return betas


betas = []
for i in range(len(perror1)):
    print(perror1[i])
    betas.append(rateofconv(perror1[i],deltax))
print("betas=",betas)





plt.figure(7)
plt.plot(-np.log(deltax), np.log(perror1))
plt.figure(7)
plt.suptitle("Rate of Convergence (Case 1)")
plt.title(r"$\alpha$=3")
plt.xlabel(r"-log($\Delta x$)")
plt.ylabel("log(error)")
plt.show()

