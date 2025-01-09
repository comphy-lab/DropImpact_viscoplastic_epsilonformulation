import numpy as np
import matplotlib.pyplot as plt

J=0.01
mu=1
mu_maxs=[1,10,100,500,1000,2000,10000]
x=np.arange(0.00001,0.002,0.00001)


plt.figure(figsize=(8,6))
for mu_max in mu_maxs:
    y=np.minimum((J/(x)+mu)*x,mu_max*x)
    plt.plot(x,y,label=rf"$muMAX={mu_max}$")
    plt.xlabel(r"||D||")
    plt.ylabel(r"$\sigma$")
    plt.legend()
    plt.savefig("model-mu.png",format="png",dpi=300)

