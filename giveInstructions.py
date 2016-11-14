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
filename = "landingInstruction.txt"
sys = AST1100SolarSystem(seed)

sys.landOnPlanet(1,filename)
