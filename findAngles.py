import numpy as np
#from PIL import Image
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
from mpl_toolkits.mplot3d import Axes3D
import seaborn



seed = 75041
sys = AST1100SolarSystem(seed)
time = 107461.01


def getPlanetAngle(time):

    rot = 1.0/(sys.period[1])*time/(3600*24.)

    if rot > 1.0:
        rot -= np.floor(rot)
    return rot*2*np.pi

def periodicAngle(angle):
    return angle%(2*np.pi)#angle - np.floor(angle/(2*np.pi))*(2*np.pi)


pos = -np.array([ 2661046.24640579, -1501372.65179696  ,      0.         ])
phi = np.arctan2(pos[1],pos[0])
theta = np.arccos(pos[2]/float(norm(pos)))
if phi < 0:
    phi += 2*np.pi

phi = periodicAngle(phi - getPlanetAngle(time))


print "r [km]: ", (norm(pos)/1000. - sys.radius[1])
print "r_tot [km]",norm(pos)/1000.
print "Theta: ",theta
print "Phi: ",phi
print "Phi in deg", np.rad2deg(phi)
