import numpy as np
import bolton
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import Subplot
from matplotlib.ticker import FuncFormatter, Formatter
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

# Variables

C_to_K = 273.15
skew_slope = 40

def parse_SPC(filename, skip_rows = 6):
    '''Reads in pressure, p, height, z, temperature, T, dewpoint, Td, and wind speed/direction from a .txt file as retrieved from SPC'''
    dtype = [('p', float), ('z', float), ('T', float), ('Td', float), ('wind_dir', float), ('wind_spd', float)]
    data = np.genfromtxt(filename, dtype = dtype, skip_header=skip_rows, delimiter=',')
    return data

def x_from_TP(T,P):
    '''Takes temperature (K) and pressure (hPa) to convert to an x coordinate from a skewed temperature and logarithmic pressure profile'''
    x = T - skew_slope*np.log(P)
    return x
def y_from_P(P):
    '''Takes pressure (hPa) to convert to an y coordinate from a logarithmic pressure profile'''
    y = -1*np.log(P)
    return y
def T_from_xP(x,P):
    '''Takes an x coordinate and pressure (hPa) from a skewed temperature and logarithmic pressure profile and returns temperature (K)'''
    T = x + skew_slope*np.log(P)
    return T
def P_from_y(y):
    '''Takes a y coordinate from a skewed temperature and logarithmic pressure profile and returns pressure (hPa)'''
    P = np.exp(-y)
    return P
def to_thermo(x,y):
    """Transform (x,y) coordinates to T in degrees Celsius and P in mb"""
    P = P_from_y(y)
    T_C = T_from_xP(x,P) - C_to_K
    return T_C, P
def from_thermo(T_C,P):
    """Transform T_C (in degrees Celisus) and P (in mb) to (x,y)"""
    y = y_from_P(P)
    x = x_from_TP(T_C+C_to_K,P)
    return x,y

# values along the bottom and left edges
P_bottom = 1050.0
P_top = 150.
T_min = -40 + C_to_K
T_max = 50 + C_to_K
x_min = x_from_TP(T_min,P_bottom)
x_max = x_from_TP(T_max,P_bottom)
y_min = y_from_P(P_bottom)
y_max = y_from_P(P_top)

# Define plotting levels for pressure, temperatures and mixing ratios
P_levels = np.arange(1000, 150-50, -50)
T_C_levels = np.arange(-80, 41, 10)
T_levels = T_C_levels + C_to_K
theta_levels = np.arange(-40+C_to_K, 100+10+C_to_K, 10)
theta_ep_levels = theta_levels.copy()
theta_e_levels = theta_levels.copy()
mixing_ratios = np.asarray([0.4, 1, 2, 3, 4, 8, 12, 16, 20])/1000.

# Convert pressure, temperatures, and mixing ratios to x and y coordinates
P_all = np.arange(P_bottom, P_top-1, -1)
y_P_levels = y_from_P(P_levels)
y_all_P = y_from_P(P_all)
x_T_levels = [x_from_TP(Ti, P_all) for Ti in T_levels]
x_thetas = [x_from_TP(bolton.theta_dry(theta_i, P_all), P_all) for theta_i in theta_levels]
x_mixing_ratios = [x_from_TP(bolton.mixing_ratio_line(P_all, mixing_ratios_i)+C_to_K, P_all) for mixing_ratios_i in mixing_ratios]

# Create a 2D field of pressure and temperature
mesh_T, mesh_P = np.meshgrid(np.arange(-60.0, T_levels.max()-C_to_K+0.1, 0.1), P_all)
theta_ep_mesh = bolton.theta_ep_field(mesh_T, mesh_P)
theta_e_mesh = bolton.equiv_potential_T(mesh_T, mesh_P)

skew_grid_helper = GridHelperCurveLinear((from_thermo, to_thermo))

# Begin the plot
fig = plt.figure(figsize=(12,12))
ax = Subplot(fig, 1, 1, 1, grid_helper = skew_grid_helper)

def format_coord(x,y):
    T, p = to_thermo(x,y)
    return "{0:5.1f} C, {1:5.1f} mb".format(float(T), float(p))
ax.format_coord = format_coord

fig.add_subplot(ax)

# Plot temperature, pressure, theta_dry, and mixing ratios for all given levels
for yi in y_P_levels:
    ax.plot((x_min, x_max), (yi,yi), color = (1., 0.8, 0.8))
for x_T in x_T_levels:
    ax.plot(x_T, y_all_P, color = (1., 0.5, 0.5))
for x_theta in x_thetas:
    ax.plot(x_theta, y_all_P, color = (1., 0.7, 0.7))
for x_mixing_ratio in x_mixing_ratios:
    good = P_all >= 600 # restrict mixing ratio lines to below 600 mb
    ax.plot(x_mixing_ratio[good], y_all_P[good], color = (0.8, 0.8, 0.6))
    
n_moist = len(theta_ep_levels)
moist_colors = ((0.6, 0.9, 0.7),)*n_moist # define colors for theta_ep contours

# Plot 2D grid of theta_ep and theta_e
ax.contour(x_from_TP(mesh_T+C_to_K, mesh_P), y_from_P(mesh_P), theta_ep_mesh, theta_ep_levels, colors = moist_colors)
ax.contour(x_from_TP(mesh_T+C_to_K, mesh_P), y_from_P(mesh_P), theta_e_mesh, theta_e_levels, colors = 'c')
ax.axis((x_min, x_max, y_min, y_max))

# Import real sounding data (as .txt)
sounding_data = parse_SPC('test_sounding.txt')
snd_T = sounding_data['T']
good_T = (snd_T > -100.) & (snd_T < 60.)
snd_P = sounding_data['p']
snd_Td = sounding_data['Td']
x_snd_T = x_from_TP(snd_T+C_to_K, snd_P)
x_snd_Td = x_from_TP(snd_Td+C_to_K, snd_P)
y_snd_P = y_from_P(snd_P)
ax.plot(x_snd_T, y_snd_P, linewidth=2, color='r')
ax.plot(x_snd_Td, y_snd_P, linewidth=2, color='g')

# Make it shine!
plt.title('LCH 181025/1200')
plt.show()

