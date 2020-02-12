#importing the necessary libraries needed for this program
import numpy as np
import astropy.units as u
from ReadFile import Read #imports necessary info from ReadFile we made

#define function that will return the total mass of any Galaxy
def ComponentMass(File, PType):

    time, total, data = Read(File)  #reading the data we need from the ReadFile

    index = np.where(data['type'] == PType)  #creates an array to store indexes of particles of desired PType
    
    mnew = data['m'][index]*1e10*u.Msun #returns mass in solar mass units of the PType

    # Mass value rounded to 3 decimal places
    Mass = np.sum(mnew)

    return Mass

# Using the above defined ComponentMass to calculate the total mass of each galaxy

# Accounting for each type of particle in order to find total mass of each galaxy using a fo rloop

# Particles in the Halo for each galaxy
MW_1 = np.round(ComponentMass('MW_000.txt',1)/1e12,3)
M31_1 = np.round(ComponentMass('M31_000.txt',1)/1e12,3)
M33_1 = np.round(ComponentMass('M33_000.txt',1)/1e12,3)

# Particles in the Disk for each galaxy
MW_2 = np.round(ComponentMass('MW_000.txt',2)/1e12,3)
M31_2 = np.round(ComponentMass('M31_000.txt',2)/1e12,3)
M33_2 = np.round(ComponentMass('M33_000.txt',2)/1e12,3)

# Particles in the Bulge for each galaxy
MW_3 = np.round(ComponentMass('MW_000.txt',3)/1e12,3)
M31_3 = np.round(ComponentMass('M31_000.txt',3)/1e12,3)
M33_3 = np.round(ComponentMass('M33_000.txt',3)/1e12,3)

# Computing total mass within each galaxy
MW_tot = MW_1 + MW_2 + MW_3
M31_tot = M31_1 + M31_2 + M31_3
M33_tot = M33_1 + M33_2 + M33_3

# Rounding the masses to three decimals
MW_tot = np.round(MW_tot,3)
M31_tot = np.round(M31_tot,3)
M33_tot = np.round(M33_tot,3)

# Computing total mass within the Local Group (LG)
LG = MW_tot + M31_tot + M33_tot

LG = np.around(LG,3) # Rounds LG to three decimals

# Calculating the Baryon functions for each galaxy and the LG
MW_stellar = MW_2 + MW_3
M31_stellar = M31_2 + M31_3
M33_stellar = M33_2 + M33_3
LG_stellar = MW_stellar + M31_stellar + M33_stellar

MW_bar = MW_stellar / MW_tot
M31_bar = M31_stellar / M31_tot
M33_bar = M33_stellar / M33_tot
LG_bar = LG_stellar / LG

# Rounding the masses to three decimals
MW_bar = np.round(MW_bar,3)
M31_bar = np.round(M31_bar,3)
M33_bar = np.round(M33_bar,3)
LG_bar = np.round(LG_bar,3)

print("MW Halo mass is ",MW_1)
print("M31 Halo mass is ",M31_1)
print("M33 Halo mass is ",M33_1)

print("MW Disk mass is ",MW_2)
print("M31 Disk mass is ",M31_2)
print("M33 Disk mass is ",M33_2)

print("MW Bulge mass is ",MW_3)
print("M31 Bulge mass is ",M31_3)
print("M33 Bulge mass is ",M33_3)

print("MW mass is ", MW_tot)
print("M31 mass is ",M31_tot)
print("M33 mass is ",M33_tot)
print("LG mass is ",LG)

print("MW Baryon function is ",MW_bar)
print("M31 Baryon function is ",M31_bar)
print("M33 Baryon function is ",M33_bar)
print("LG Baryon function is ",LG_bar)
