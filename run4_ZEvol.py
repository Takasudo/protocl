from math import *
import numpy as np
import sys
args = sys.argv

# [0] setting

snap = [21, 25, 33, 50, 67, 98]
redshift = [4.0, 3.0, 2.0, 1.0, 0.5, 0.01]
age = [1.540, 2.145, 3.285, 5.878, 8.587, 13.667] # [Gyr]  
calclabel = 'tmp'

# mass bin

z0idlist_bin1 = [i for i in range (0,34)] # mbin=0 : mass > 14.5
#z0idlist_bin1 = [i for i in range (35,98)] # mbin=1 : 14.25 < mass < 14.5
#z0idlist_bin1 = [i for i in range (98,252)] # mbin=2 : 14.0 < mass < 14.25

mbin=0

if mbin==0:
    mnote='       M > $10^{14.5} M_\odot$\n        (N='+str(len(z0idlist_bin1))+')'
elif mbin==1:
    mnote='$10^{14.25} < $M < $10^{14.5} M_\odot$\n (N='+str(len(z0idlist_bin1))+')'
elif mbin==2:
    mnote='$10^{14} < $M < $10^{14.25} M_\odot$\n (N='+str(len(z0idlist_bin1))+')'

print(z0idlist_bin1) 

# [1] load data

N = len(snap)

co2ph = np.zeros(N)

R_core_bin1 = np.zeros((len(z0idlist_bin1),N))
R_cl_bin1 = np.zeros((len(z0idlist_bin1),N))
R_100_bin1 = np.zeros((len(z0idlist_bin1),N))
SFR_total_bin1 = np.zeros((len(z0idlist_bin1),N))  # SFR (all)
SFR_core_bin1 = np.zeros((len(z0idlist_bin1),N))   # SFR (R < R_core)
f_SFR_core_bin1 = np.zeros((len(z0idlist_bin1),N)) # Fractional SFR (R < R_core)

i = 0
for z0id in z0idlist_bin1:
    j = 0
    for snapnum in snap:
        npz = np.load('npzfiles/proto_radial_data_FixRmax_snap' + str(snap[j]) + "_z0id" + str(z0idlist_bin1[i]) + "_" + str(calclabel) + ".npz")
        co2ph[j] = npz['conv']
        R_core_bin1[i][j] = npz['size_core']
        R_cl_bin1[i][j] = npz['size_50']
        R_100_bin1[i][j] = npz['size_100']
        #SFR_total_bin1[i][j] = max(npz['sfr_proto'][npz['radius_list']<R_cl_bin1[i][j]])
        SFR_total_bin1[i][j] = max(npz['sfr_proto'][npz['radius_list'] < R_100_bin1[i][j]])
        SFR_core_bin1[i][j] = max(npz['sfr_proto'][npz['radius_list']<R_core_bin1[i][j]])
        if SFR_total_bin1[i][j]>0:
            f_SFR_core_bin1[i][j] = SFR_core_bin1[i][j] / SFR_total_bin1[i][j]
        j +=1
    i += 1

# [2] calc properties

R_core_mean_bin1 = np.zeros(N)
R_cl_mean_bin1 = np.zeros(N)
SFR_core_mean_bin1 = np.zeros(N)
SFR_total_mean_bin1 = np.zeros(N)
f_SFR_mean_bin1 = np.zeros(N)

for j in range (N):
    R_core_mean_bin1[j] += np.median(R_core_bin1[:,j])
    R_cl_mean_bin1[j] += np.median(R_cl_bin1[:,j])
    #SFR_core_mean_bin1[j] += np.median(SFR_core_bin1[:,j])
    #SFR_total_mean_bin1[j] +=np.median(SFR_total_bin1[:,j])
    SFR_core_mean_bin1[j] += np.mean(SFR_core_bin1[:,j])
    SFR_total_mean_bin1[j] +=np.mean(SFR_total_bin1[:,j])   
    f_SFR_mean_bin1[j] += np.median(f_SFR_core_bin1[:,j])

# [3] fitting function

def Rproto_fit(t_gyr, a, b, c):
    return a * (1 + t_gyr/b)**(-c)

def Rcore_fit(t_gyr, a, b):
    return a * (1 + t_gyr)**b

def fcore_fit(t_gyr, a, b, c):
    return a * (1 + t_gyr/b)**c

def SFRproto_fit(t_gyr, a, b, c):
    return a * (1 + t_gyr/b)**(-c)

def SFRcore_fit(t_gyr, a, b, c):
    return a * (1 + t_gyr/b)**(-c)

# fit param (c : core, p : proto)

# bin 0, 50% radius
'''
Rc1_bin1, Rc2_bin1 = 0.3e3, 0.6
Rp1_bin1, Rp2_bin1, Rp3_bin1 = 1.5e4,10.0,3.0

Sc1_bin1, Sc2_bin1, Sc3_bin1 = 0.5e3,4.0,3.0
Sp1_bin1, Sp2_bin1, Sp3_bin1 = 2.e4,9.0,8.0
'''

# bin 0, 80% radius  
'''
Rc1_bin1, Rc2_bin1 = 0.3e3, 0.6
Rp1_bin1, Rp2_bin1, Rp3_bin1 = 1.5e4,12.0,2.5

Sc1_bin1, Sc2_bin1, Sc3_bin1 = 0.5e3,4.0,3.0
Sp1_bin1, Sp2_bin1, Sp3_bin1 = 3.e4,5.0,6.0
'''

# bin 0, 100% radius 

Rc1_bin1, Rc2_bin1 = 0.3e3, 0.6
Rp1_bin1, Rp2_bin1, Rp3_bin1 = 1.5e4,10.0,3.0

Sc1_bin1, Sc2_bin1, Sc3_bin1 = 1.e3,4.0,3.0
Sp1_bin1, Sp2_bin1, Sp3_bin1 = 2.e4,8.0,6.0

# bin 1
'''
Rc1_bin1, Rc2_bin1 = 0.24e3, 0.6
Rp1_bin1, Rp2_bin1, Rp3_bin1 = 1.2e4,10.0,3.0

Sc1_bin1, Sc2_bin1, Sc3_bin1 = 0.5e3,4.0,3.0
Sp1_bin1, Sp2_bin1, Sp3_bin1 = 1.e4,6.0,5.5
'''

# bin2
'''
Rc1_bin1, Rc2_bin1 = 0.2e3, 0.6
Rp1_bin1, Rp2_bin1, Rp3_bin1 = 1.e4,10.0,3.0

Sc1_bin1, Sc2_bin1, Sc3_bin1 = 0.5e3, 4.0, 3.0
Sp1_bin1, Sp2_bin1, Sp3_bin1 = 0.9e4, 7.0, 7.0
'''

# [4] Cosmic SFR Density
'''
import astropy.units as u
from astropy.cosmology import Planck13, z_at_value

def psi_from_t(t_gyr):
    # [M_sun / yr / Mpc^3]
    z = z_at_value(Planck13.age, t_gyr * u.Gyr)
    return 0.015 * (1+z)**2.7 / (1 + ((1+z)/2.9)**5.6)
psi_from_t_vec = np.vectorize(psi_from_t)
timelist = np.linspace(0.1,13.7,100)
sfrd = psi_from_t_vec(timelist)
'''

# [5] plot

import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

fig = plt.figure(figsize=(5,7))

ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

for i in range(len(z0idlist_bin1)):
    ax1.plot(np.array(age), R_core_bin1[i], marker='x', lw=0, ms=5, color='red', alpha=0.3)
    ax1.plot(np.array(age), R_cl_bin1[i], marker='+', lw=0, ms=5, color='blue', alpha=0.3)
    ax2.plot(np.array(age), SFR_total_bin1[i], marker='x',lw=0, ms=5, color='blue', alpha=0.3)
    ax2.plot(np.array(age), SFR_core_bin1[i], marker='+', lw=0, ms=5, color='red', alpha=0.3)

ax1.plot(np.array(age), R_core_mean_bin1, marker='x', markeredgewidth=2.5, lw=0, ms=9, color='red', label='Core Radius')
ax1.plot(np.array(age), R_cl_mean_bin1, marker='+', markeredgewidth=2.5, lw=0, ms=9, color='blue', label='Proto Cluster Radius')
ax2.plot(np.array(age), SFR_core_mean_bin1, marker='x', markeredgewidth=2.5, lw=0, ms=9, color='red', label='SFR (<$R_{core}$')
ax2.plot(np.array(age), SFR_total_mean_bin1, marker='+', markeredgewidth=2.5, lw=0, ms=9, color='blue', label='SFR (all)')

ax1.plot(np.array(age),Rcore_fit(np.array(age[:]),Rc1_bin1, Rc2_bin1), lw=2, color='red', ls=':', alpha=0.6)
ax1.plot(np.array(age),Rproto_fit(np.array(age[:]),Rp1_bin1, Rp2_bin1, Rp3_bin1), lw=2, color='blue', ls=':', alpha=0.6)
ax2.plot(np.array(age),SFRcore_fit(np.array(age[:]),Sc1_bin1, Sc2_bin1, Sc3_bin1), lw=2, color='red', ls=':', alpha=0.6)
ax2.plot(np.array(age),SFRproto_fit(np.array(age[:]),Sp1_bin1, Sp2_bin1, Sp3_bin1 ), lw=2, color='blue', ls=':', alpha=0.6)

# compare with R_50, mbin=0
'''
ax1.plot(np.array(age),Rcore_fit(np.array(age[:]),0.3e3, 0.6), lw=1, color='red', ls=':', alpha=0.4)
ax1.plot(np.array(age),Rproto_fit(np.array(age[:]),1.5e4,10.0,3.0), lw=1, color='blue', ls=':', alpha=0.4)
ax2.plot(np.array(age),SFRcore_fit(np.array(age[:]),0.5e3,4.0,3.0), lw=1, color='red', ls=':', alpha=0.4)
ax2.plot(np.array(age),SFRproto_fit(np.array(age[:]),2.e4,9.0,8.0), lw=1, color='blue', ls=':', alpha=0.4)
'''

# Cosmic SFRD
'''
ax2.plot(timelist, sfrd*1e4, lw=1, color='gray', ls=':', alpha=0.75)
ax2.text(10,5e2,'Cosmic\nSFR',color='gray', alpha=0.75)
'''

ax1.set_xlabel("cosmic age [Gyr]")
ax1.set_ylabel("radius [ckpc/h]")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.legend()

ax2.set_xlabel("cosmic age [Gyr]")
ax2.set_ylabel("SFR [M$_\odot$ yr$^{-1}$]")
ax2.set_ylim(3,3e4) 
ax2.set_xscale("log")
ax2.set_yscale("log")

ax2.text(4,5e3,mnote,fontsize=12.5)

ax1.set_xlim(1,17)
ax2.set_xlim(1,17)

plt.savefig("fig_z_protocr_"+str(calclabel)+"_m"+str(mbin)+".pdf",bbox_inches='tight')

# [5] SFR density

from scipy.special import erf, hyp2f1

def rho0f0(r,t,s_0):
    # t in [Gyr]
    # r in [ckpc/h]
    SFRcore = SFRcore_fit(t,Sc1_bin1, Sc2_bin1, Sc3_bin1)
    Rcore = Rcore_fit(t,Rc1_bin1, Rc2_bin1)
    #
    ## (1) : exp(-r**2)
    #tmp1 = SFRcore / (pi * s_0**3.0 * Rcore**3.0)
    #tmp2 = sqrt(pi) * erf(1/s_0) - 2.0 * exp(-1/s_0**2.0)/s_0
    #radial = exp(-(r/(s_0*Rcore))**2.0)
    #result = tmp1/tmp2 * radial
    #
    ## (2) : exp(-r)
    tmp1 = SFRcore / (4.0 * pi * s_0**3.0 * Rcore**3.0)
    tmp2 = 2.0 - ((1/s_0)*(1/s_0 + 2.0) + 2.0) * exp(-1/s_0)
    radial = exp(-r/(s_0*Rcore))
    result = tmp1 / tmp2 * radial
    #
    return result

def rho1f1(r,t,s_1):
    # t in [Gyr]  
    # r in [ckpc/h]
    SFRcore = SFRcore_fit(t,Sc1_bin1, Sc2_bin1, Sc3_bin1)
    SFRproto = SFRproto_fit(t,Sp1_bin1, Sp2_bin1, Sp3_bin1)
    Rcore =  Rcore_fit(t,Rc1_bin1, Rc2_bin1)
    Rproto = Rproto_fit(t,Rp1_bin1, Rp2_bin1, Rp3_bin1)
    #
    Rcut = s_1 * Rcore
    #Rcut = Rcore + (Rproto - Rcore) / s_1
    #
    ## (1) : exp(-r**2)
    ## result 1 : normalization from r = 0 to Rcut
    ## result 2 : normalization from r = 0 to infty
    #r_tmp = Rproto / Rcut
    #tmp1 = (SFRproto - SFRcore)/ (pi * Rcut**3.0)
    #tmp2 = sqrt(pi) * erf(r_tmp) - 2.0 * exp(-r_tmp**2.0)*r_tmp
    #radial = exp(-(r/Rcut)**2.0)
    #result1 = tmp1/tmp1 * radial
    #result2 = tmp1 / 4.0 / 0.443 * radial
    #
    ## (2) : exp(-r)
    tmp1 = (SFRproto - SFRcore)/ (4.0 * pi * Rcut**3.0)
    radial = exp(-(r / Rcut))
    result = tmp1 / 2.0 * radial
    #
    return result

rho0vec = np.vectorize(rho0f0)
rho1vec = np.vectorize(rho1f1)

radius_list = np.logspace(0,4.5,50)
rs = radius_list**3
