import numpy as np

def FDM(n, L, alpha, case, order):
    deltax = L/n
    if case == 1:
        K = np.zeros((n-1,n-1))
        f = np.zeros((n-1)) 
        if order == 2:
            kappa = 2+alpha**2*deltax**2
        elif order == 4:
            beta = alpha**2*deltax**2/12
            kappa = (2+10*beta)/(1-beta)
        elif order == 6:
            kappa = 2+alpha**2*deltax**2 + alpha**4*deltax**4/12 + alpha**6*deltax**6/360
        elif order == 8:
            kappa = 2+alpha**2*deltax**2 + alpha**4*deltax**4/12 + alpha**6*deltax**6/360 + alpha**8*deltax**8/20160
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
    elif case == 2:
        K = np.zeros((n,n))
        f = np.zeros((n))
        if order == 2:
            kappa = 2+alpha**2*deltax**2
        elif order == 4:
            beta = alpha**2*deltax**2/12
            kappa = (2+10*beta)/(1-beta)
        elif order == 6:
            kappa = 2+alpha**2*deltax**2 + alpha**4*deltax**4/12 + alpha**6*deltax**6/360
        elif order == 8:
            kappa = 2+alpha**2*deltax**2 + alpha**4*deltax**4/12 + alpha**6*deltax**6/360 + alpha**8*deltax**8/20160
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