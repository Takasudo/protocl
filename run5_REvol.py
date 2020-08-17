
from run4_ZEvol import *

from math import *
import numpy as np
from scipy.stats import binned_statistic
import sys
args = sys.argv

snap = [21, 25, 33, 50, 67, 98]
redshift = [4.0, 3.0, 2.0, 1.0, 0.5, 0.01]
age = [1.540, 2.145, 3.285, 5.878, 8.587, 13.667] # [Gyr]
calclabel = 'tmp'

# mass bin

z0idlist = [i for i in range (0,34)]
mbin = 0

import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

N = len(snap)
Rbins = 25
cl_radius = np.array([[np.zeros(Rbins) for k1 in range (len(z0idlist))] for l1 in range (N)])
cl_sfrd = np.array([[np.zeros(Rbins) for k2 in range (len(z0idlist))] for l2 in range (N)])

#--- load data

j=0
for z0id in z0idlist:
    i = 0
    for snapnum in snap:
        #npz = np.load('npzfiles/proto_radial_data_snap'+str(snapnum)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
        npz = np.load('npzfiles/proto_radial_data_FixRmax_snap' + str(snapnum) + "_z0id" + str(z0id) + "_" + str(calclabel) + ".npz")
        radius = (npz['radius_list'][1:] + npz['radius_list'][:-1]) / 2.0
        dr = np.diff(npz['radius_list'])
        sfr_shell = np.diff(npz['sfr_proto'])
        sfr_density = sfr_shell / (4.0 * pi * radius ** 2.0 * dr)
        #
        sfrd_bin, bin_edge, bin_num = binned_statistic(np.log(radius),sfr_density, statistic='mean', bins=Rbins)
        dlnr = bin_edge[1] - bin_edge[0]
        bin_lnr = bin_edge[1:] - dlnr/2.
        bin_r = np.exp(bin_lnr)
        #
        cl_radius[i][j] += bin_r
        cl_sfrd[i][j] += sfrd_bin
        #
        i+=1
    j += 1

#--- stack

cl_sfrd_mean = np.array([np.zeros(Rbins) for l3 in range (N)])

for j in range (len(z0idlist)):
    for i in range (len(snap)):
        cl_sfrd_mean[i] += cl_sfrd[i][j] / float(len(z0idlist))

#--- plot

fig, ax = plt.subplots(N,1,figsize=(8,20))

list_score = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
list_scl = [5, 5, 5, 5, 5, 5]

j=0
for i in range(N):
    j = 0
    for z0id in z0idlist:
        if j<10:
            cl = 'gray'
            al = 0.9 * (1+j)**(-0.5)
        elif j<20:
            cl = 'gray'
            al = 0.9 * (1+(j-10))**(-0.5)
        else:
            cl= 'gray'
            al = 0.8 * (1+(j-20))**(-0.5)
        al = 0.3

        ax[i].set_xscale("log")
        ax[i].set_yscale("log")
        ax[i].set_xlim(1e1, 3e4)
        ax[i].set_ylim(1e-3, 1e3)
        
        ax[i].plot(cl_radius[i][j], cl_radius[i][j] ** 3.0 * cl_sfrd[i][j], c=cl, lw=0, marker='+', ms=5, alpha=al)
        ax[i].plot(cl_radius[i][j], cl_radius[i][j] ** 3.0 * cl_sfrd_mean[i], lw=0, c='black', marker='+',ms=15, markeredgewidth=2, alpha=0.5)

        # for legend
        ax[i].plot(cl_radius[i][j - 1], cl_radius[i][j - 1] ** 3.0 * cl_sfrd_mean[i], lw=0, c='black', marker='+', ms=15,
                 markeredgewidth=2, label='Illustris')

        # rho.pdf from z_evol.py

        tgyr = age[i]

        r_core = Rcore_fit(tgyr, Rc1_bin1, Rc2_bin1)
        r_cl = Rproto_fit(tgyr, Rp1_bin1, Rp2_bin1, Rp3_bin1)

        ax[i].plot(radius_list, rs * rho0vec(radius_list, tgyr, list_score[i]), ls='-.', color='red', label='core',lw=1.5)
        ax[i].plot(radius_list, rs * rho1vec(radius_list, tgyr, list_scl[i]), ls='-.', color='blue', label='CL', lw=1.5)

        j+=1

ax[i].set_xlabel("Radius r [ckpc/h]",fontsize=13)
ax[i].set_ylabel(r"$r^3 \dot{\rho}_*$ [M$_\odot$ yr$^{-1}$]", fontsize=13)
#ax[i].legend()

plt.savefig("fig_r_protocr_"+str(calclabel)+"_m"+str(mbin)+".pdf",bbox_inches='tight')
