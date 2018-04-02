#! /usr/bin/python3
from numpy import *
import matplotlib.pyplot as plt
import seaborn as sns

# --------------------------
#             ^
#             |
#
#         d; V={4,7,10}; grid(1000)
#
#             |
#             ^
# --------------------------
# delta, phi = incline; distort
# k11, k22, k33 = 12.6E-12, 6.1E-12, 18.65E-12
# eps_para, eps_vert = 12.23*eps0, 4.75*eps0
# the_0 = 5 * (pi / 180)
# phi_0 = 270 * (pi / 180)
# ddp = 0.75


#%%
eps0 = 8.854E-12

# k11, k22, k33 = 12.6E-12, 6.1E-12, 18.65E-12
# eps_para, eps_vert = 12.23*eps0, 4.75*eps0

k11, k22, k33 = 13.7E-12, 7.0E-12, 16.8E-12
eps_para, eps_vert = 7.3*eps0, 3.6*eps0

deps = eps_para - eps_vert

the_0 = 5 * (pi / 180)
phi_0 = 270 * (pi / 180)
h = 1E-12 # 1um
NUM_GRID = 1000
d = h*NUM_GRID

p_0 = d / 0.75
q_0 = 2*pi / p_0

#%%
# ---- init ------
V_spec = [4, 7, 10]
V_spec = [10]

THE = 0.5*pi*random.rand(NUM_GRID)
THE[0] = the_0
PHI = 0.5*pi*phi_0*random.rand(NUM_GRID)
PHI[0] = phi_0


def f(the):
    return k11*cos(the)**2 + k33*sin(the)**2

def g(the):
    return (k22*cos(the)**2 + k33*sin(the)**2)*cos(the)**2

def df(the, the_p, the_n):
    return (k33 -k11)*sin(2*the)*(the_p - the_n) / (2*h)

def dg(the, the_p, the_n):
    return (-k22 + (-k22 + k33)*cos(2*the))*sin(2*the)*(the_p - the_n)/(2*h)

#%%
# --- iter ----
# for V in V_spec:
#     THE_RES, PHI_RES = [], []
#     for i in range(1, NUM_GRID-1):
#         # THE, PHI = new_THE, new_PHI
#         the, the_p, the_n = THE[i],THE[i+1],THE[i-1]
#         phi, phi_p, phi_n = PHI[i],PHI[i+1],PHI[i-1]
#         while True: # -- iterate n
#             the_src, phi_src = the, phi
#
#             the = (0.125 / f(the)) * ( 0 + 4*f(the)*(the_p + the_n) - 2*h*k22*q_0*sin(2*the)*(phi_p - phi_n) + 2*h*(the_p - the_n)*df(the, the_p ,the_n))
#
#             phi = (0.25 / the)*(h*k22*q_0*sin(2*the)*(the_p - the_n) + 2*g(the)*(phi_p + phi_n) + h*dg(the, the_p ,the_n)*(phi_p - phi_n))
#
#             if abs(the - the_src) < 1E-9 and abs(phi - phi_src) < 1E-9:
#                 THE_RES.append(the)
#                 PHI_RES.append(phi)
#                 break

for V in V_spec:
    while True: # -- iterate n
        THE_SRC = THE
        PHI_SRC = PHI
        for i in range(1, NUM_GRID-1):
            # THE, PHI = new_THE, new_PHI
            the, the_p, the_n = THE_SRC[i],THE_SRC[i+1],THE_SRC[i-1]
            phi, phi_p, phi_n = PHI_SRC[i],PHI_SRC[i+1],PHI_SRC[i-1]

            the = (0.125 / f(the)) * ( 0 + 4*f(the)*(the_p + the_n) - 2*h*k22*q_0*sin(2*the)*(phi_p - phi_n) + 2*h*(the_p - the_n)*df(the, the_p ,the_n))

            phi = (0.25 / the)*(h*k22*q_0*sin(2*the)*(the_p - the_n) + 2*g(the)*(phi_p + phi_n) + h*dg(the, the_p ,the_n)*(phi_p - phi_n))

            THE[i], PHI[i] = the, phi
        if max(abs(THE - THE_SRC)) < 1E-9 and max(abs(PHI - PHI_SRC)) < 1E-9:
            break

    # THE[1:-1] = THE_RES
    # PHI[1:-1] = PHI_RES
    plt.figure()
    plt.plot(THE)
    plt.figure()
    plt.plot(PHI)
    plt.show()
