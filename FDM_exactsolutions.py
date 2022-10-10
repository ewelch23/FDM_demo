import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt

## exact solutions
## for a bar of length L, consider two boundary cases:
## 1; U(0)=0, U(L)=0 the temperature at the left end of the bar is 0, and the temperature at the right end is 100 C.
## 2; U'(0)=0, U(L)=100 the bulk heat flux at the left end is 0, and the temperature at the right end is 100 C.
## The exact solutions for temperature and heat flux for both these cases are plotted below.
L = 1
A = 1
alphas = [0.2,0.4,0.75,3,5,7]
u_exact1 = lambda x,alpha: 100/np.sinh(alpha * L) * np.sinh(alpha * x)
u_exact2 = lambda x,alpha: 100/np.cosh(alpha * L) * np.cosh(alpha * x)
du_exact1 = lambda x,alpha: (100*alpha*np.cosh(alpha*x))/np.sinh(L*alpha)
du_exact2 = lambda x,alpha: (100*alpha*np.sinh(alpha*x))/np.cosh(L*alpha)
Q_exact1 = lambda x,alpha: A*du_exact1(x,alpha)
Q_exact2 = lambda x,alpha: A*du_exact2(x,alpha)
x_exactplot = np.linspace(0,L)

for i in range(len(alphas)):
    plt.figure(1)
    plt.plot(x_exactplot, u_exact1(x_exactplot,alphas[i]),label=r"$\alpha$="+str(alphas[i]))
    plt.figure(2)
    plt.plot(x_exactplot, u_exact2(x_exactplot,alphas[i]),label=r"$\alpha$="+str(alphas[i]))
    plt.figure(3)
    plt.plot(x_exactplot,Q_exact1(x_exactplot,alphas[i]),label=r'$\alpha$='+str(alphas[i]))
    plt.figure(4)
    plt.plot(x_exactplot,Q_exact2(x_exactplot,alphas[i]),label=r'$\alpha$='+str(alphas[i]))

plt.figure(1)
plt.title("Exact Analytical Solution: Temperature, Case 1")
plt.xlabel("x")
plt.ylabel("Temperature (C)")
plt.legend()
plt.savefig("hw1_exactTemp_case1.png")

plt.figure(2)
plt.title("Exact Analytical Solution: Temperature, Case 2")
plt.xlabel("x")
plt.ylabel("Temperature (C)")
plt.legend()
plt.savefig("hw1_exactTemp_case2.png")

plt.figure(3)
plt.title('Exact Analytical Solution: Heat Flux, Case 1')
plt.xlabel('x')
plt.ylabel('Heat Flux (Q)')
plt.legend()
plt.savefig('hw1_exactQ_case1.png')

plt.figure(4)
plt.title('Exact Analytical Solution: Heat Flux, Case 2')
plt.xlabel('x')
plt.ylabel('Heat Flux (Q)')
plt.legend()
plt.savefig('hw1_exactQ_case2.png')

#################################################################