#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

def tar(x):
    return x - np.cos(x)

def div2solve(func, bd, eps=1E-5):
    k = 0
    a, b = bd
    X = []
    while True:
        k = k+ 1
        x = (a + b) / 2
        X.append(x)
        if np.sqrt(func(x)**2) < eps: break
        if func(a)*func(x) < 0:
            b = x
        elif func(a)*func(x) > 0:
            a = x
    return (x, X)

bd = (-10, 10)
plt.figure()
rg = np.linspace(bd[0], bd[1], 1000)
plt.plot(rg, tar(rg))

#%%
plt.figure()
(x, X) = div2solve(tar, bd)
plt.plot(X)
plt.title(x)
plt.show()
