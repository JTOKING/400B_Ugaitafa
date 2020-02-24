# Homework 5
# Mass Distribution of Each Galaxy
# Justin Ugaitafa

# import modules
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterOfMass import CenterOfMass
from astropy.constants import G # imports gravitaional constant

# Class to define the mass distribution of a given galaxy and snapshot
class MassProfile:

    def __init__(self, File, PType):
    # Initialize the instance of this class with the following properties:

        # reads data in given file using Read
        self.time, self.total, self.data = Read(File)

        #creates an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == PType)

        # add a string of the filenumber to the value "000"
        ilbl = '000' + str(Snap)
        # remove all but last 3 digits
        ilbl = ilbl[-3:]
        self.File = "%s_"%(galaxy) + ilbl + '.txt'
    
        # store the mass, positions, velocities of only the particles of the
        # given type into arrays
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]*u.kpc
        self.y = self.data['y'][self.index]*u.kpc
        self.z = self.data['z'][self.index]*u.kpc

        self.gname = galaxy

    def MassEnclosed(self, PType, r):
    # Function that will compute the mass within any given radius with respect
    # to the COM for a specified galaxy and specified component of that galaxy.
    # input: the array of all the radii of the total mass enclosed
    # return: the array of the mass of each particle within each specified
    #         radius.
        # read in COM position from last HW
        COM = CenterOfMass(COM_P)
        Position = COM.COM_P(0.1)
        # 1) find the x,y,z positions of each indicidual particle within
        #    specified radius.
        x = self.x - Position
        y = self.y - Position
        z = self.z - Position
        # 2) find the mass of each particle
        m = self.m
        # computes magnitude of position vector
        r = np.sqrt(x**2 + y**2 + z**2)

        # for loop to iterate over the radius given at each array element
        for i in range(0,r):
            particlepos = np.where(r1 <= r[i]) # individual particle position
            particlemass = m[particlepos] # individual particle mass
            mass[i] = sum(particlemass) # total mass of particles

        return mass*1e10*u.Msun
    
    def MassEnclosedTotal(self, r):
    # Function to compute mass enclosed within an enclosed radius
    # input: array(x,y,z) of radii
    # return: array of masses (Msun) representing total mass(bulge+disk+halo) at
    #         each radius of the input array
        EncHaloM = self.MassEnclosed(0,r) # Mass of particles enclsoed in halo
        EncDiskM = self.MassEnclosed(1,r) # Mass of particles enclosed in disk
        EncBulgeM = self.MassEnclosed(2,r) # Mass of particles enclosed in bulge

        totalM = EncHaloM + EncDiskM + EncBulgeM
        return totalM

    def HernquistMass(self, r, a, Mhalo):
    # Function that returns Hernquist mass profile
    # inputs:
        # r = distance from galaxy (kpc)
        # a = scale radius (kpc)
        # Mhalo = total dark matter halo mass (10^12 Msun)
    # return:
        # total dark matter enclosed within r (Msun)
        return np.around(Mhalo*r**2/(a+r)**2,2)*1e12*u.Msun

    def CircularVelocity(self, PType, r):
    # Function that finds the circular speeds in (km/s), rounded to 2 decimal
    # places
    # inputs:
        # PType = particle type
        # array with radii
    # return:
        # array of circular speeds
        mass = self.MassEnclosed(PType,r) # mass of each particle type
        CircVel = np.sqrt(G*mass / r*u.kpc) # circ velocity of each particle
        return np.around(CircVel,2)

    def CircularVelocityTotal(self, r):
    # Function that finds the total circular velocity created by each galaxy
    # component at each radius of the radii array
    # input:
        # array of radii
    # return:
        # array of circular velocity
        M = self.MassEnclosedTotal(r) # total particle mass within each radius
        CircVel2 = np.sqrt(G*M / r*u.kpc) # total circ vel of the particles
        return np.around(CircVel2,2)

    def HernquistVCirc(self, r, a, Mhalo):
    # Function that finds circular speed (km/s) within the Hernquist M profile
    # inputs:
        # r = distance from galaxy (kpc)
        # a = scale radius (kpc)
        # Mhalo = total dark matter halo mass (10^12 Msun)
    # return:
        # circular speed (km/s) of mass in Hernquist M profile
        M2 = self.HernquistMass(r, a, Mhalo)
        CircVel3 = np.sqrt(G*M2 / r*u.kpc)
        return np.around(CircVel3,2)

# Defining the mass profile for each galaxy
MWMP = MassProfile('MW',0)
M31MP = MassProfile('M31',0)
M33MP = MassProfile('M33',0)
R = np.arrange(0.1,30,0.5)

# Defining the masses of each galaxy component and that galaxy's total mass
MWHalo = MWMP.MassEnclosed(0,R)
MWDisk = MWMP.MassEnclosed(1,R)
MWBulge = MWMP.MassEnclosed(2,R)
MWTotal = MWMP.MassEnclosedTotal(R)

M31Halo = M31MP.MassEnclosed(0,R)
M31Disk = M31MP.MassEnclosed(1,R)
M31Bulge = M31MP.MassEnclosed(2,R)
M31Total = M31MP.MassEnclosedTotal(R)

M33Halo = M33MP.MassEnclosed(0,R)
M33Disk = M33MP.MassEnclosed(1,R)
M33Bulge = M33MP.MassEnclosed(2,R)
M33Total = M33MP.MassEnclosedTotal(R)

# There's a call error where I defined the mass profiles but I have no idea how
# to fix it.  That's why I have no graphs.
