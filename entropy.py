#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
#%%
# k = 0
nx, nt = 1000, 2000

kTv = 0.01; # temp

eps = 1E-3

P = 0.99

T_spec = np.linspace(0 + 0.01, 0 + 0.22 , nt)
P_res = []
for kTv in T_spec:
    P_iter = []
    # print(kTv)
    while True:
        dx = 1/nx # each slice
        # print(dx)
        nume, deli = 0, 0

        for k in range(nx):
            x = k*dx # current loc
            ak = 1.5*x**2 - 0.5
            bk = ak / (kTv)
            nume += ak*np.exp(bk*P)
            deli += np.exp(bk*P)
        P_src = P
        P = nume / deli

        P_iter.append(P)
        if np.abs(P_src - P) < eps:
            P_res.append(P)
            break;
    # plt.figure()
    # plt.plot(P_iter)
    # plt.show()
#%%
for k in np.arange(0,len(T_spec)):
    print("% 5.5f"%T_spec[k], " : ", "% 5.5f"%P_res[k])
#%%
plt.figure()
plt.plot(T_spec, P_res)
plt.show()
