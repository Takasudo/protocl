import illustris_python as il
import sys
from math import *
import numpy as np

args = sys.argv

snap = 99
basePath = '/Volumes/ExtremeSSD/WorkingFolder/proto/TNG300-1/output'

#subhalos = il.groupcat.loadSubhalos(basePath,snap,fields=['SubhaloGrNr','SubhaloPos'])
FoFgroup = il.groupcat.loadHalos(basePath,snap,fields=['GroupFirstSub','GroupPos','GroupNsubs','GroupMass'])

# Generate sh file

outputfilename = args[1]

file = open(outputfilename+".sh",'w')
snap = ["21 4 ", "25 3 ","33 2 ","50 1 ","67 0.5 ","98 0.01 ","99 0.0 "]

mass, i, j = 15, 0, 0
mass_min, mass_max = 14.5, 15.5

file.write("# note : \n")
file.write("# mass_min = " + str(mass_min)+"\n")
file.write("# mass_max = " + str(mass_max)+"\n")
file.write("# \n")

while (mass > mass_min):
    mass = np.log10(FoFgroup['GroupMass'][i]) + 10
    if (mass >= mass_min) and (mass <= mass_max):
    #if i==102 or i==103:  # one needs some manual modification because FoFgroup['GroupMass'][i] is not exactly ordered by mass
        for z in snap:
            file.write("python3 run2_ProtoCL.py "+ z + str(i) + " "+ str(j) + " " + str(j+FoFgroup['GroupNsubs'][i]-1) + "\n")
            file.write("python3 run3_RadialSubhalo.py " + z + str(i) + "\n")

    j +=FoFgroup['GroupNsubs'][i]
    i +=1

file.write("echo 'edit and run run4_ZEvol.py'")
file.close()
