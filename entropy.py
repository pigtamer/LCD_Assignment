#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
#%%
# k = 0
nx, nt = 500, 500

eps = 1E-3

liP = [0.99, -0.99, 0.01]

RESes = []
T_spec = [  np.linspace(0 + 0.01, 0 + 0.28 , nt),
            np.linspace(0 + 0.01, 0 + 0.28 , nt),
            np.linspace(0 + 0.01, 0 + 0.28 , nt)
            ]
for kP in range(len(liP)):
    P = liP[kP]
    P_res = []
    for kTv in T_spec[kP]:
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
            if np.abs(P_src - P) < eps:
                P_res.append(P)
                break
    RESes.append(P_res)
    print("Init point: % 5.3f"%P)
    for k in np.arange(0,len(T_spec[kP])):
        print("% 5.5f"%T_spec[kP][k], " : ", "% 5.5f"%P_res[k])

#%%
plt.figure(figsize = (10,7))
plt.plot(T_spec[0], RESes[0], '-', T_spec[1], RESes[1], '-',T_spec[2], RESes[2], '-')
plt.show()
