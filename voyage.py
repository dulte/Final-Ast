import numpy as np
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
import seaborn
import sys


seed = 75041
sys = AST1100SolarSystem(seed)

sys.sendSatellite("instructions.txt")
