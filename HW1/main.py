# ---------------------------------#
# -- main script for NSE 565 HW1 --#
# --------- Austin Warren ---------#
# ---------- Winter 2022 ----------#
# ---------------------------------#

import numpy as np
import matplotlib.pyplot as plt

# define central difference scheme function to use for each case
def cds_ss(num_volumes, tot_length, velocity, density, diffusion, source, left, right):
    """Function to perform central difference scheme to solve one-dimensional steady state transport with convection and diffusion.

    Parameters
    ----------
    num_volumes : float
        The number of discretized volumes.
    tot_length : float
        The total length of the pipe in meters.
    velocity : float
        The average velocity of the flow in meters per second.
    density : float
        The density of the flow in kilograms per cubic meter.
    diffusion : float
        The diffusion coefficient in kilogram-seconds per meter.
    source : float
        The source term in [units].
    left : float
        The left boundary condition.
    right : float
        The right boundary condition.

    Returns
    -------
    phi : numpy.ndarray
        Solved flux profile.
    """
    dx = tot_length/num_volumes
    phi = np.zeros(num_volumes)
    A = np.zeros((num_volumes,num_volumes))
    Q = np.zeros(num_volumes)

    

    for j in range(num_volumes):

        if j == 0:
            A[j,0] = density*velocity/2 + 3*diffusion/dx
            A[j,1] = density*velocity/2 - diffusion/dx
            Q[j] = density*velocity*left + 2*diffusion*left/dx

        elif j == num_volumes-1:
            A[j,j-1] = -density*velocity/2 - diffusion/dx
            A[j,j] = -density*velocity/2 + 3*diffusion/dx
            Q[j] = -density*velocity*right + 2*diffusion*right/dx

        else:
            A[j,j-1] = -density*velocity/2 - diffusion/dx
            A[j,j] = 2*diffusion/dx
            A[j,j+1] = density*velocity/2 - diffusion/dx
            Q[j] = 0

    phi = np.linalg.solve(A,Q)
    return phi



# define analytical solution function to use for each case
def analytical(num_volumes, tot_length, velocity, density, diffusion, left, right):
    """Function to analytically solve one-dimensional steady state transport with convection and diffusion.

    Parameters
    ----------
    num_volumes : float
        The number of discretized volumes.
    tot_length : float
        The total length of the pipe in meters.
    velocity : float
        The average velocity of the flow in meters per second.
    density : float
        The density of the flow in kilograms per cubic meter.
    diffusion : float
        The diffusion coefficient in kilogram-seconds per meter.
    left : float
        The left boundary condition.
    right : float
        The right boundary condition.

    Returns
    -------
    phi : numpy.ndarray
        Solved flux profile.
    """
    x = np.linspace(0, tot_length, num=num_volumes,endpoint=True)
    phi = left + (right-left) * (np.exp(density*velocity*x/diffusion)-1)/(np.exp(density*velocity*tot_length/diffusion)-1)
    return phi



#
def error_calc(phi_exact, phi):
    """Function to calculate error.

    Parameters
    ----------
    phi_exact : numpy.ndarray
        The exact solution.
    phi : numpy.ndarray
        The approximate solution.

    Returns
    -------
    error : float
        The average error for the approximate solution.
    """
    error = np.sum(np.abs(phi_exact - phi)) / len(phi)
    return error


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


error = np.zeros(len(velocity))
# solve
for k in range(len(velocity)):
    phi = cds_ss(num_volumes[k], tot_length, velocity[k], density, diffusion, source, left, right)
    phi_exact = analytical(num_volumes[k], tot_length, velocity[k], density, diffusion, left, right)
    error[k] = error_calc(phi_exact,phi)
    
    
    x = np.linspace(0, tot_length, num=num_volumes[k])

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    plt.figure(facecolor='w', edgecolor='k', dpi=200)
    plt.plot(x, phi, '-k', label='CDS')
    plt.plot(x, phi_exact, '-r', label='Analytical')
    plt.figlegend(loc='right', bbox_to_anchor=(0.4,0.2))
    plt.grid(b=True, which='major', axis='both')
    plt.savefig('HW1/plots/graph_case'+str(k+1)+'.png',transparent=True)
    plt.savefig('HW1/plots/graph_case'+str(k+1)+'.svg',transparent=True)



# generate latex table
out_file = open('HW1/tabs/error_tab.tex','w')
out_file.write(
                '\\begin{table}[htbp]\n'+
                '\t \centering\n'+
                '\t \caption{Error values for each case.}\n'+
                '\t \\begin{tabular}{cc}\n'+
                '\t\t \\toprule\n'+
                '\t\t Case & Error \\\ \n'+
                '\t\t \midrule \n'+
                '\t\t 1 & '+str(error[0])+' \\\ \n'+
                '\t\t 2 & '+str(error[1])+' \\\ \n'+
                '\t\t 3 & '+str(error[2])+' \\\ \n'+
                '\t\t \\bottomrule \n'+
                '\t \end{tabular} \n'+
                '\t \label{tab:error} \n'+
                '\end{table}'
)


