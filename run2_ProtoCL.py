# --- create_subhalo_npz.py :
# (1) read protocluster subhalo data from illustris simulation
# (2) define protocluster membership radius (d_N)
# (3) save all data to "npzdata_protocluster_subhalo" file

import sys
from math import *
import numpy as np
import illustris_python as il

args = sys.argv

flag_gal = 0 # [0: from snapshot, 1: from npz file]

calclabel = 'tmp'

snap = int(args[1])
redshift = float(args[2])
z0id, sub_start, sub_end = int(args[3]), int(args[4]), int(args[5])

sim_box_size = 205000.0
basePath = '/Volumes/ExtremeSSD/WorkingFolder/proto/TNG300-1/output'

# Cosmology

hubble = 0.7

# (1) set subhalo data

# (1a) load data

subhalo_pos = []
subhalo_mass = []
subhalo_sfr = []
subhalo_size = []
subhalo_bh = []
subhalo_id = []
subhalo_gr = []
subhalo_grRc200 = []
subhalo_grRc500 = []
subhalo_grRm200 = []
subhalo_grRh200 = []
subhalo_grPos = []

for i in range (sub_start, sub_end):
    tree = il.sublink.loadTree(basePath,99,i,fields=['SubfindID','SubhaloGrNr','SubhaloPos','SnapNum','SubhaloMass','SubhaloSFR','SubhaloHalfmassRad','SubhaloBHMdot','Group_R_Crit200','Group_R_Crit500','Group_R_Mean200','Group_R_TopHat200','GroupPos'],onlyMPB=False)
    #--- select subhalos with tree['SnapNum']==33 (z=2)"
    if tree == None:
        continue
    else:
        for n in range(len(tree['SnapNum'])):
            snapnum = tree['SnapNum'][n]
            if(snapnum==snap):
                subhalo_pos.append(tree['SubhaloPos'][n])
                subhalo_sfr.append(tree['SubhaloSFR'][n])
                subhalo_mass.append(tree['SubhaloMass'][n])
                subhalo_size.append(tree['SubhaloHalfmassRad'][n])
                subhalo_bh.append(tree['SubhaloBHMdot'][n])
                subhalo_id.append(tree['SubfindID'][n])
                subhalo_gr.append(tree['SubhaloGrNr'][n])
                subhalo_grRc200.append(tree['Group_R_Crit200'][n])
                subhalo_grRc500.append(tree['Group_R_Crit500'][n])
                subhalo_grRm200.append(tree['Group_R_Mean200'][n])
                subhalo_grRh200.append(tree['Group_R_TopHat200'][n])
                subhalo_grPos.append(tree['GroupPos'][n])
    del tree


np.savez('npzfiles/subhalo_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel),
         subhalo_pos = subhalo_pos,
         subhalo_sfr = subhalo_sfr,
         subhalo_mass = subhalo_mass,
         subhalo_size = subhalo_size,
         subhalo_bh = subhalo_bh)

# (1b) set np arrays

Nhalo = len(subhalo_mass)     # total number of subhalos
mass = np.array(subhalo_mass) # array of subhalo mass
sfr = np.array(subhalo_sfr)   # array of subhalo SFR
pos_x = np.array([row[0] for row in subhalo_pos]) # x of subhalos
pos_y = np.array([row[1] for row in subhalo_pos]) # y of subhalos
pos_z = np.array([row[2] for row in subhalo_pos]) # z of subhalos

print("Sub halo num : ",Nhalo)

# (2) properties of subhalos within different radii

# (2a) find central halo

# most massive subhalo

for i in range (len(mass)):
    if mass[i] == np.max(mass): 
        most_massive_grPos = subhalo_grPos[i]
        most_massive_grRc200 = subhalo_grRc200[i]
        most_massive_grRc500 = subhalo_grRc500[i]
        most_massive_grRm200 = subhalo_grRm200[i]
        most_massive_grRh200 = subhalo_grRh200[i]

# (2b) define center and core

center = most_massive_grPos
coreR = most_massive_grRh200

# (2c) define radius  (d20, d50, d80, d100)   

diff = np.absolute(subhalo_pos - np.array([center for n in range(Nhalo)])) # List of |(subhalo position) - (center)|
distance_from_center = []

Lcrit = sim_box_size * 0.7
diff[diff>Lcrit] = sim_box_size - diff[diff>Lcrit]

for d in range (len(diff)):
    dist = (diff[d][0]**2.0 + diff[d][1]**2.0 + diff[d][2]**2.0)**0.5
    distance_from_center.append(dist)
distance = np.array(distance_from_center) # List of distances between subhalo_pos and center

sorted_distance = np.sort(distance)       # Sorted distances
d_20 = sorted_distance[int(Nhalo*0.2)]    # 20% distance
d_50 = sorted_distance[int(Nhalo*0.5)]    # 50% distance
d_80 = sorted_distance[int(Nhalo*0.8)]    # 80% distance
d_100 = sorted_distance[Nhalo-1]          # 100% distance

if len(subhalo_size)==1:
    d_100 = subhalo_size[0]  # re-define 100% distance

xmin, xmax = center[0] - d_100, center[0] + d_100
ymin, ymax = center[1] - d_100, center[1] + d_100
zmin, zmax = center[2] - d_100, center[2] + d_100 


# (2d) define physical quantities within d20, d50, d80, d100

subhalo_d20_pos = []
subhalo_d20_sfr = []
subhalo_d20_size = []
subhalo_d20_mass = []

subhalo_d50_pos = []
subhalo_d50_sfr = []
subhalo_d50_size = [] 
subhalo_d50_mass = []

subhalo_d80_pos = []
subhalo_d80_sfr = []
subhalo_d80_size = []
subhalo_d80_mass = []

subhalo_d100_pos = []
subhalo_d100_sfr = []
subhalo_d100_size = []
subhalo_d100_mass = []

for i in range (Nhalo):
    if (distance[i] < d_20):
        subhalo_d20_pos.append(subhalo_pos[i])
        subhalo_d20_sfr.append(subhalo_sfr[i])
        subhalo_d20_size.append(subhalo_size[i])
        subhalo_d20_mass.append(subhalo_mass[i])
    if (distance[i] < d_50):
        subhalo_d50_pos.append(subhalo_pos[i])
        subhalo_d50_sfr.append(subhalo_sfr[i])
        subhalo_d50_size.append(subhalo_size[i])
        subhalo_d50_mass.append(subhalo_mass[i])
    if (distance[i] < d_80):
        subhalo_d80_pos.append(subhalo_pos[i])
        subhalo_d80_sfr.append(subhalo_sfr[i])
        subhalo_d80_size.append(subhalo_size[i])
        subhalo_d80_mass.append(subhalo_mass[i])
    if (distance[i] < d_100):
        subhalo_d100_pos.append(subhalo_pos[i])
        subhalo_d100_sfr.append(subhalo_sfr[i])
        subhalo_d100_size.append(subhalo_size[i])
        subhalo_d100_mass.append(subhalo_mass[i])

print("Load subhalo data")

# (2e) from comiving kpc/h to physical Mpc

d_20_pMpc = d_20*1.e-3/(1.+redshift)/hubble    # from [ckpc/h] to [Mpc]
d_50_pMpc = d_50*1.e-3/(1.+redshift)/hubble
d_80_pMpc = d_80*1.e-3/(1.+redshift)/hubble
d_100_pMpc = d_100*1.e-3/(1.+redshift)/hubble

print(" from ckpc/h to kpc : ",1./(1.+redshift)/hubble)

# (3) save data

np.savez("npzfiles/protocluster_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel),
         CoM = center,
         xmin = xmin,
         ymin = ymin,
         zmin = zmin,
         xmax = xmax,
         ymax = ymax,
         zmax = zmax,
         d_20 = d_20,
         d_50 = d_50,
         d_80 = d_80,
         d_100 = d_100,
         d20_physMpc = d_20_pMpc,
         d50_physMpc = d_50_pMpc,
         d80_physMpc = d_80_pMpc,
         d100_physMpc = d_100_pMpc,
         subhalo_d20_sfr = subhalo_d20_sfr, # [M_sun/yr]
         subhalo_d50_sfr = subhalo_d50_sfr,
         subhalo_d80_sfr = subhalo_d80_sfr,
         subhalo_d100_sfr = subhalo_d100_sfr,
         subhalo_d20_mass = subhalo_d20_mass, # [10^10 M_sun/h]
         subhalo_d50_mass = subhalo_d50_mass,
         subhalo_d80_mass = subhalo_d80_mass,
         subhalo_d100_mass = subhalo_d100_mass,
         subhalo_d20_pos = subhalo_d20_pos,
         subhalo_d50_pos = subhalo_d50_pos,
         subhalo_d80_pos = subhalo_d80_pos,
         subhalo_d100_pos = subhalo_d100_pos,
         core_size = coreR
         )

''' 

Previous Note

# Snap and redshift

 99 : z=0
 33 : z=2
 50 : z=1
 67 : z=0.5
 91 : z=0.1
 98 : z=0.01

# local fof info

#z0id = 0       # put FoF Halo ID
#sub_start = 0  # put subhalo start
#sub_end = 1907 # put subhalo end

#z0id, sub_start, sub_end = 0, 0, 1907
#z0id, sub_start, sub_end = 1, 2066, 3201
#z0id, sub_start, sub_end = 2, 3202, 4089
#z0id, sub_start, sub_end = 3, 4090, 4975
#z0id, sub_start, sub_end = 4, 4976, 5894
#z0id, sub_start, sub_end = 5, 5895, 6626
#z0id, sub_start, sub_end = 6, 6628, 7094

properites of z=0 FoF halo

0 2066 15.291694164276123
1 1136 15.0163893699646
2 888 14.92448902130127
3 886 14.894059181213379
4 919 14.893155574798584
5 732 14.859557151794434
6 863 14.857169151306152
7 847 14.805735111236572
8 735 14.766603469848633
9 635 14.739508152008057
10 497 14.710947036743164
11 620 14.694820880889893
12 565 14.694688320159912
13 655 14.690198421478271
14 357 14.629459381103516
15 571 14.632087707519531
16 470 14.629456996917725
'''

