# ---------------------------------#
# -- main script for NSE 565 HW1 --#
# --------- Austin Warren ---------#
# ---------- Winter 2022 ----------#
# ---------------------------------#

import numpy as np
import matplotlib.pyplot as plt

# define central difference scheme function to use for each case
def cds_ss(dx, tot_length, velocity, density, diffusion, source, left, right):



# define analytical solution function to use for each case
def analytical(tot_length, velocity, density, diffusion, left, right):
    




# define constants
tot_length = 1  # pipe length in meters
density = 1     # density in kg/m^3
diffusion = 0.1 # diffusion coefficient in kg-s/m
source = 0      # source
left = 100      # left hand boundary condition
right = 50      # right hand boundary condition

# define variables
velocity = np.array([0.1, 2.5, 2.5])
num_volumes = np.array([5,5,20])
dx = tot_length / num_volumes