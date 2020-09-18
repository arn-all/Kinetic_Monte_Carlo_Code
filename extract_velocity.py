import pyvista as pv
import glob
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import yaml
import natsort
import seaborn as sns


def files_list(pattern):
    """Gives a natsorted list of files that match a pattern."""
    return natsort.natsorted(glob.glob(pattern))

def parse_config(file):
    """Parse KMC code input files. Returns a dict of key-values pairs."""
    with open(file, 'r') as f:
        return yaml.safe_load(f.read().replace("=",": ").replace("\t", " "))

def get_time(file):
    """find the current time in a vtk file (in seconds)"""
    with open(file,'r') as f:
        line = f.readlines()[1]
    return float(re.findall('\(t=(.+)\)$', line)[0])

def average_line_position(img):
    """Compute the average position along x of a dislocation, using a weighted average of the position of each segment"""
    prev_pt = img[0]
    segments = []
    position = 0
    
    def pos(num):
        # If "prev_pt" is after "point", ignore the segment
        if num<0:
            return 0
        else :
            return num

    for i, point in enumerate(img):
        if i == 0:
            # skip first iteration
            continue
        if np.isclose(point[0], prev_pt[0]):
            # if the points have the same x position
            position += point[0] * pos(point[2]-prev_pt[2])
        prev_pt = point
    return position/max(img[:,2]) # the sum divided by the cell width

def collect_positions(flist):
    """Returns the positions of all dislocations of pointdefects of a simulation."""
    collect_elems = []
    times = []

    for file in flist:
        update_mesh = pv.read(file)
        if "dislocations" in file:
            collect_elems.append(average_line_position(np.asarray(update_mesh.points)))
        elif "pointdefects" in file :
            collect_elems.append(np.asarray(update_mesh.points))
        else:
            raise RuntimeError("{} is neither dislocations or pointdefects.")
        times.append(get_time(file))
    
    return np.asarray(collect_elems), times


def speed_schmid(params):
    """Analytical expression of dislocation velocity. Non schmid expression from Zhao et al. 2018"""
    h = np.sqrt(6)/3 * float(params['lattice_param']) # size of a valley. (m)
    b = np.sqrt(3)/2 * float(params['lattice_param']) # (m)
    w0 = float(params['kpwidth_w0']) #2.31 # (burgers)
    
    a1 = float(params["nonschmid_a1"]) # 0 #0.938
    a2 = float(params["nonschmid_a2"]) #0 #0.71
    a3 = float(params["nonschmid_a3"]) #0 #4.43
    nonschmid_tc = float(params["nonschmid_tc"]) # 2.03
    
    deltaH0 = float(params['kpenergy_deltaH0']) # 1.63 # (eV)
    k = 8.617e-5 # eV.K^-1
    T = float(params["temperature"])
    stress = float(params['stress_yz'])
    
    nu = float(params['attempt_frequency']) #  9.1e11 # (s^-1)

    # Adjustable params (no dimension)
    chi = 0
    m = 0.5
    n = 0.15
    c = 2.02
    p = 0.86
    q = 1.69
    non_glide = float(params['stress_yy'])-float(params['stress_xx'])
    
    sigma1 = stress * np.cos(chi)
    sigma2 = stress * np.cos(chi + np.pi/3)
    sigma3 = non_glide * np.sin(2*chi)
    sigma4 = non_glide * np.sin(2*chi + 2*np.pi/3.)
    
    def f(x):
        return 2./(1+np.exp(2*x))
    
    
    def s_(stress, chi):
        return (sigma1+a1*sigma2) / (nonschmid_tc*float(params['peierls_stress'])*f((a2*sigma3 + a3*sigma4)/(nonschmid_tc*float(params['peierls_stress']))))
    
    def kink_pair_sep(s):
        # denoted w(s)
        return w0*(s**(-m)+c)*(1-s)**(-n)

    def enthalpy(s):
        # DeltaH(tau)
        return deltaH0*(1-s**p)**q

    S = s_(stress, chi)
    v = h * nu* (params["line_length"]-kink_pair_sep(S))/1. * np.exp(-enthalpy(S)/(k*T))
    return v/b


