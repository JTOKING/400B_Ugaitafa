# Homework 6
# Seperation and Relative Velocities of the COM of Simulated Galaxies
# Justin Ugaitafa

# import modules
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib

# importing from files from previous HWs
from ReadFile import Read
from CenterOfMass import CenterOfMass
from astropy.constants import G # imports gravitational constant

class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy and simulation snapshot
    
    
    def __init__(self, File, PType):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(File)                                                                                             

        #creates an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == PType)

        # store the mass, positions, velocities of only the particles of the
        # given type
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)

        # xcomponent Center of mass
        Acom = np.sum(a*m) / np.sum(m)
        # ycomponent Center of mass
        Bcom = np.sum(b*m) / np.sum(m)
        # zcomponent Center of mass
        Ccom = np.sum(c*m) / np.sum(m)
        
        return Acom, Bcom, Ccom
    
    
    def COM_P(self, delta, VolDec):
    # Function to specifically return the center of mass position and velocity      # inputs:                                                                       #        particle type (1,2,3)                                                  #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                    #       3D coordinates of the center of mass position (kpc)                                                                                                 

        # my first guess at the COM position by calling COMdefine
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        
        # compute the magnitude of the COM position vector.
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)


        # iterative process to determine the center of mass change reference frame to COM frame                                                                         # compute the difference between particle coordinates and the first guess at COM position
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

        # find the max 3D distance of all particles from the guessed COM                # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/VolDec
        
        # My initial value for the change in COM position between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                  # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW <= RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                         # compute the center of mass position using the particles in the reduced radius
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            
            # computes the new 3D COM position
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position and the new one.
            CHANGE = np.abs(RCOM - RCOM2)
            
            # uncomment the following line if you wnat to check this                        #print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM reduces the volume by a factor of 2 again
            RMAX = RMAX/VolDec
            
            # check this by uncommenting the following
            #print ("RMAX = ", RMAX)

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # creates a vector to store the COM position
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units using astropy and round all values and then return the COM positon vector
        COMP = np.around(COMP,2)*u.kpc
        return COMP

    def COM_V(self, XCOM,YCOM,ZCOM):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns: 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determines the position of all particles relative to the center of mass position
        xV = self.x - XCOM.value
        yV = self.y - YCOM.value
        zV = self.z - ZCOM.value
        RV = np.sqrt(xV**2 + yV**2 + zV**2)

        # determine the index for those particles within the max radius
        indexV = np.where(RV <= RVMAX.value)

        # determine the velocity and mass of those particles within the mas radius
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew =  self.m[indexV]
        
        # computes the center of mass velocity using those particles
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew,mnew)

        # create a vector to store the COM velocity
        # set the correct units using astropy and round all values to two places
        COMV = [VXCOM, VYCOM, VZCOM]
        COMV2 = np.around(COMV,2)*(u.km/u.s)

        # return the COM vector                                                                                        
        return COMV2
