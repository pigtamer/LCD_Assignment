#! /usr/bin/python3
from numpy import *
from numpy.random import *
import matplotlib.pyplot as plt
import seaborn as sns
import time

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

h = 1E-6 # 1um
NUM_GRID = 100
d = h*NUM_GRID

p_0 = d / 0.75
q_0 = 2*pi / p_0

#%%
# ---- init randomly ------
the_0, phi_0, vol_0 = 5 * (pi / 180), 0 * (pi / 180), 4

THE = the_0*rand(NUM_GRID)
PHI = phi_0*rand(NUM_GRID)
VOL = vol_0*rand(NUM_GRID)
THE[0], PHI[0], VOL[0] = the_0, phi_0, vol_0
THE[-1] = the_0
PHI[-1] = 90*(pi/180)
# THE = the_0 * ones(NUM_GRID)
# PHI = phi_0 * ones(NUM_GRID)
# VOL = vol_0 * ones(NUM_GRID)

THE_NXT, PHI_NXT, VOL_NXT  = THE.copy(), PHI.copy(), VOL.copy()

#%%
# ---- Assistant Func ----
def f(the):
    return k11*cos(the)**2 + k33*sin(the)**2

def g(the):
    return (k22*cos(the)**2 + k33*sin(the)**2)*cos(the)**2

def df(the, the_p, the_n):
    return (k33 -k11)*sin(2*the)*(the_p - the_n) / (2*h)

def dg(the, the_p, the_n):
    return (-k22 + (-k22 + k33)*cos(2*the))*sin(2*the)*(the_p - the_n)/(2*h)

def calc_theta(the, the_p, the_n, vol_p, vol_n, phi_p, phi_n):
    return((0.125 / f(the)) * ( deps*cos(the)*sin(the)*(vol_p - vol_n)**2 + 4*f(the)*(the_p + the_n) - 2*h*k22*q_0*sin(2*the)*(phi_p - phi_n) + 2*h*(the_p - the_n)*df(the, the_p ,the_n)))

def calc_phi(the, the_p, the_n, vol_p, vol_n, phi_p, phi_n):
    return((0.25 / g(the))*(h*k22*q_0*sin(2*the)*(the_p - the_n) + 2*g(the)*(phi_p + phi_n) + h*dg(the, the_p ,the_n)*(phi_p - phi_n)))

def calc_volt(the, the_p, the_n, vol_p, vol_n, phi_p, phi_n):
    return((2*(eps_para * sin(the)**2 + eps_vert * cos(the)**2 )*(vol_p + vol_n) + deps*cos(the)*sin(the)*(the_p - the_n)*(vol_p - vol_n)) / (4*(eps_para * sin(the)**2 + eps_vert * cos(the)**2)))



#%%
# --- iter ----
tic = time.time()
thres = 1E-6
idx = 0
while True: # -- iterate n
    idx += 1
    # if idx > 1:
    THE, PHI, VOL = THE_NXT.copy(), PHI_NXT.copy(), VOL_NXT.copy()
    for k in range(1, NUM_GRID-1): #update
        # XX_p, XX_n : positive, negative. XX(k+1), XX(k-1)
        the, the_p, the_n = THE[k],THE[k+1],THE[k-1]
        phi, phi_p, phi_n = PHI[k],PHI[k+1],PHI[k-1]
        vol, vol_p, vol_n = VOL[k],VOL[k+1],VOL[k-1]

        THE_NXT[k] = calc_theta(the, the_p, the_n, vol_p, vol_n, phi_p, phi_n)

        PHI_NXT[k] = calc_phi(the, the_p, the_n, vol_p, vol_n, phi_p, phi_n)

        VOL_NXT[k] = calc_volt(the, the_p, the_n, vol_p, vol_n, phi_p, phi_n)

    if sum(abs(THE - THE_NXT)) < thres and sum(abs(PHI - PHI_NXT)) < thres and sum(abs(VOL - VOL_NXT)) < thres: break
    # if idx > 5000: break
toc = time.time() - tic

print(">> Summary <<\n| %5.3f s ; %4d rds |\n Init Phi:\t%5.3f, %5.3f (rad)\n      Delta:\t%5.3f, %5.3f (rad)\n"%(toc, idx, PHI[0], PHI[-1], THE[0], THE[-1]))


plt.figure()
plt.plot(THE_NXT*(180/pi))
plt.title("Delta")
plt.figure()
plt.plot(PHI_NXT*(180/pi))
plt.title("Phi")
plt.figure()
plt.plot(VOL_NXT)
plt.title("Voltage")
plt.show()
