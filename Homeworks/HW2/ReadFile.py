import numpy as np
import astropy.units as u

def Read(File): #defining the Read file as a function

    file = open(File,'r')

    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr

    line2 = file.readline()
    label, value = line2.split()
    particle_num = float(value)

    file.close()
    
    data = np.genfromtxt(File,dtype=None,names=True,skip_header=3)

    return time,particle_num,data
  
time,particle_num,data = Read('MW_000.txt')
#time,particle_num,data = Read('M31_000.txt')
#time,particle_num,data = Read('M33_000.txt')

#print(data['m'][1])
