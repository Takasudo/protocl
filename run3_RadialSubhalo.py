# --- radial.py 

from math import *
import numpy as np
import matplotlib.pyplot as plt
import sys
args = sys.argv

snap = int(args[1])
redshift = float(args[2])

z0id = int(args[3])
calclabel = 'tmp' 

hubble = 0.7

sim_box_size = 205000.0

note = 'FixRmax'

# (1) load npz files

npz_subhalo = np.load("npzfiles/protocluster_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
npz_halos = np.load('npzfiles/subhalo_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")

# (2) prepare for unit conversion

length_kpc = 1./(1.+redshift)/hubble # from [ckpc/h] to [kpc]
km_to_cm = 1.e5
kpc_to_cm = 3.086e21

# (3a) subhalo : center and size of protocluster

center = npz_subhalo['CoM']

d_20 = npz_subhalo['d_20']
d_50 = npz_subhalo['d_50']
d_80 = npz_subhalo['d_80']
d_100 = npz_subhalo['d_100']

print(d_20, d_50, d_80, d_100, sim_box_size)

# (3b) subhalo : position and properties

subhalo_pos = npz_halos['subhalo_pos']
subhalo_sfr = npz_halos['subhalo_sfr']
subhalo_mass = npz_halos['subhalo_mass']
subhalo_size = npz_halos['subhalo_size']
subhalo_bh = npz_halos['subhalo_bh']

subhalo_dx = np.absolute(np.array(subhalo_pos[:,0] - center[0]))
subhalo_dy = np.absolute(np.array(subhalo_pos[:,1] - center[1]))
subhalo_dz = np.absolute(np.array(subhalo_pos[:,2] - center[2]))

Lcrit = sim_box_size * 0.7

subhalo_dx[subhalo_dx>Lcrit] = sim_box_size - subhalo_dx[subhalo_dx>Lcrit]
subhalo_dy[subhalo_dy>Lcrit] = sim_box_size - subhalo_dy[subhalo_dy>Lcrit]
subhalo_dz[subhalo_dz>Lcrit] = sim_box_size - subhalo_dz[subhalo_dz>Lcrit]

subhalo_distance_from_center = (subhalo_dx**2.0 + subhalo_dy**2.0 + subhalo_dz**2.0)**0.5

SFhalo_distance = subhalo_distance_from_center[subhalo_sfr>0]
highSFhalo_loc = subhalo_distance_from_center[subhalo_sfr>10]
highSFhalo_sfr = subhalo_sfr[subhalo_sfr>10]
subhalo_mass_sum = np.sum(subhalo_mass)
core_radius = npz_subhalo['core_size']

print("load subhalo (N = ",len(subhalo_mass),")")

# (3c) all illustris subhalo

#'''
import illustris_python as il

basePath = '/Users/taka/Desktop/work_proto/TNG300-2/output/'
allsubhalos = il.groupcat.loadSubhalos(basePath,snap,fields=['SubhaloPos','SubhaloSFR','SubhaloMass','SubhaloBHMdot'])
allsubhalo_pos = allsubhalos['SubhaloPos']
allsubhalo_sfr = allsubhalos['SubhaloSFR']
allsubhalo_mass = allsubhalos['SubhaloMass']
allsubhalo_bh = allsubhalos['SubhaloBHMdot']
allsubhalo_distance_from_center = (np.array(allsubhalo_pos[:,0] - center[0])**2.0 + np.array(allsubhalo_pos[:,1] - center[1])**2.0 + np.array(allsubhalo_pos[:,2] - center[2])**2.0)**0.5

print("load all subhalo in snapshot")
#'''

# (4) properties within radius R

#'''
R_bins = 500
r_min = 1.0
#r_max = 1.e4 + max(subhalo_distance_from_center)
r_max = 5.e4
radius_list = np.logspace(log10(r_min),log10(r_max),R_bins) # [ckpc/h]
volume_within_R = 4.0 * pi * (radius_list*length_kpc*kpc_to_cm)**3.0 / 3.0 # [cm^3]

sum_of_subhalo_sfr_within_R = np.zeros(R_bins)
sum_of_subhalo_mass_within_R = np.zeros(R_bins)
sum_of_allsubhalo_sfr_within_R = np.zeros(R_bins)
sum_of_allsubhalo_mass_within_R = np.zeros(R_bins)
sum_of_allsubhalo_bh_within_R = np.zeros(R_bins)
sum_of_subhalo_bh_within_R = np.zeros(R_bins)
num_of_subhalo_within_R = np.zeros(R_bins)
num_of_sf_subhalo_within_R = np.zeros(R_bins)

for r in range(R_bins):
    R = radius_list[r]

    # value of X for each subhalo within R

    subhalo_sfr_within_R = subhalo_sfr[(subhalo_distance_from_center<R)]
    subhalo_mass_within_R = subhalo_mass[(subhalo_distance_from_center<R)]
    subhalo_sf_within_R = subhalo_sfr[(subhalo_distance_from_center<R)&(subhalo_sfr>0)]
    allsubhalo_sfr_within_R = allsubhalo_sfr[(allsubhalo_distance_from_center<R)]
    allsubhalo_mass_within_R = allsubhalo_mass[(allsubhalo_distance_from_center<R)]
    allsubhalo_bh_within_R = allsubhalo_bh[(allsubhalo_distance_from_center<R)]
    subhalo_bh_within_R = subhalo_bh[(subhalo_distance_from_center<R)]

    # take sum of all subhalos within R

    sum_of_subhalo_sfr_within_R[r] = sum(subhalo_sfr_within_R)
    sum_of_subhalo_mass_within_R[r] = sum(subhalo_mass_within_R) 
    sum_of_subhalo_bh_within_R[r] = sum(subhalo_bh_within_R)
    sum_of_allsubhalo_sfr_within_R[r] = sum(allsubhalo_sfr_within_R)
    sum_of_allsubhalo_mass_within_R[r] = sum(allsubhalo_mass_within_R)
    sum_of_allsubhalo_bh_within_R[r] = sum(allsubhalo_bh_within_R)
    sum_of_subhalo_bh_within_R[r] = sum(subhalo_bh_within_R)
    num_of_subhalo_within_R[r] = len(subhalo_sfr_within_R)
    num_of_sf_subhalo_within_R[r] = len(subhalo_sf_within_R)
#'''

# (5) take diff

'''
shell_width = np.diff(radius_list)
mean_radius = (radius_list[1:] + radius_list[0:-1])/2.0

total_in_shell_volume = np.diff(volume_within_R)                            # V : volume of spherical shell
total_in_shell_sim_particle_num = np.diff(num_of_sim_particle_within_R)     # N : gas particle num within the shell

sum_in_shell_gas_proton_num = np.diff(sum_of_gas_proton_num_within_R)
sum_in_shell_gas_mag_energy = np.diff(sum_of_gas_mag_energy_within_R)
sum_in_shell_gas_turbB = np.diff(sum_of_gas_turbB_energy_within_R)
sum_in_shell_gas_sfr = np.diff(sum_of_gas_sfr_within_R)
sum_in_shell_gas_volume = np.diff(sum_of_gas_volume_within_R) # Note : This should not be different from shell volume (total_in_shell_volume)

sum_in_shell_dm_mass = np.diff(sum_of_dm_mass_within_R)
sum_in_shell_SFhalo_num = np.diff(sum_of_SFhalo_within_R)

'''

# (6) plot and save

#'''

plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

plt.plot(radius_list, sum_of_subhalo_sfr_within_R,c='gray',label='SFR (protocluster subhalo)',lw=4)
plt.plot(radius_list, sum_of_allsubhalo_sfr_within_R,c='cyan',label='SFR (all subhalo)',lw=2) 
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Radius r [ckpc/h]", fontsize=15)
plt.ylabel("SFR (<r)", fontsize=15)
plt.legend(fontsize=13)

plt.title('log$_{10}$(M)='+str(np.log10(subhalo_mass_sum*1e10)),fontsize=12)

plt.plot([d_50,d_50],[1e2,1e3],c='blue',ls=':',alpha=0.7,lw=2)
plt.plot([core_radius,core_radius],[1e2,1e3],c='red',ls='--',alpha=0.7)

plt.xlim(1e2,5e4)
plt.savefig('tmpfig/tmpfig_proto_radialdata_'+note+'_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+'.pdf')

np.savez('npzfiles/proto_radial_data_'+note+'_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel),
         radius_list = radius_list, # [ckpc/h]
         r_max = r_max,
         conv = length_kpc, # from [ckpc/h] to [kpc]
         z = redshift,
         sfr_proto = sum_of_subhalo_sfr_within_R,
         sfr_all = sum_of_allsubhalo_sfr_within_R,
         mass_proto = sum_of_subhalo_mass_within_R,
         mass_all = sum_of_allsubhalo_mass_within_R,
         bh_proto = sum_of_subhalo_bh_within_R,
         bh_all = sum_of_allsubhalo_bh_within_R,
         size_core = core_radius,
         size_50 = d_50,
         size_80 = d_80,
         size_100 = d_100)

#'''
