import numpy as np
import matplotlib.pyplot as plt

J=0.01
mu=1
epsilons=[1,0.1,0.01,0.001,0.0001,0.00001,0.000001]
x=np.arange(0.00001,0.002,0.00001)
n=1

plt.figure(figsize=(8,6))
for epsilon in epsilons:
    y=(J/(x+epsilon)+mu*pow(x+epsilon,n-1))*x
    plt.plot(x,y,label=rf"$\epsilon={epsilon}$")
    plt.xlabel(r"||D||")
    plt.ylabel(r"$\sigma$")
    plt.legend()
    plt.savefig("model-e.png",format="png",dpi=300)

