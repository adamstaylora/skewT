import numpy as np


# constants

C_to_K = 273.15
c_P_dry = 1005.7
c_v_dry = 718
eps = 0.6220
k_dry = 0.2854
c_L = 4181.3
L_v = 2.5E-6

def sat_vapor_pressure(T):
    '''Takes a temperature (C) and returns saturation vapor pressure (hPa) via eqn 10 (Bolton, 1980)'''
    e_s = 6.112*np.exp((17.67*T)/(T+243.5))
    return e_s
def sat_vapor_temperature(e_s):
    '''Takes a saturation vapor pressure (hPa) and returns a temperature (C) via eqn 11 (Bolton, 1980)'''
    T = (243.5*np.log(e_s)-440.8)/(19.48-np.log(e_s))
    return T
def sat_mixing_ratio(P, T):
    '''Takes a pressure (hPa) and temperature (C) and returns a mixing ratio (kg/kg) via 
    w = epsilon*e_s/(P-e_s) from AMS Glossary of Meteorology Definition'''
    w = eps*sat_vapor_pressure(T)/(P-sat_vapor_pressure(T)) 
    return w
def mixing_ratio_line(P, w_s):
    '''Takes a pressure (hPa) and mixing ratio (kg/kg) and returns temperature (C) via mixing ratio definition and eqn 11 (Bolton, 1980)'''
    e_s = (w_s*P)/(w_s+eps) 
    T = sat_vapor_temperature(e_s) 
    return T
def RH(T, P, w):
    '''Takes a pressure (hPa), temperature (C), and mixing ratio (kg/kg) and returns relative humidity (%) '''
    RH = 100*w/sat_mixing_ratio(P, T) 
    return RH
def T_LCL(T, RH):
    '''Takes a temperature (C) and relative humidity (%) and returns the temperature (K) at the lifting condensation level via eqn 22 (Bolton, 1980)'''
    A = 1/(T+ C_to_K-55)
    B = np.log(RH/100.0)/2840.0
    T = 1/(A-B)
    return T
def theta_dry(theta, P, P_0 = 1000.):
    '''Takes a potential temperature (K) and returns the dry potential temperature (K) via eqn 23 (Bolton, 1980)'''
    T = theta*(P_0/P)**(-1*k_dry) 
    return T
def pseudoeq_potential_T(T, P, w, P_0 = 1000.0):
    '''Takes a temperature (C), pressure (hPa), mixing ratio (kg/kg), and reference pressure (hPa) and returns the pseudoequivalent potential temperature (K) via eqn 43 (Bolton, 1980)'''
    rh = RH(T, P, w) 
    T_lcl = T_LCL(T, rh)
    A = ((3.376/T_lcl)-0.00254)*w*1000*(1+0.81*(w))
    theta_ep = np.exp(A)*(T+C_to_K)*(P_0/P)**(0.2854*(1-0.28*(w)))
    return theta_ep
def theta_ep_field(T, P, P_0 = 1000.0):
    '''Takes temperature (C), pressure (hPa), and reference pressure (hPa) and returns pseudoequivalent potential temperature (K), eliminating the need to input mixing ratio'''
    w = sat_mixing_ratio(P,T)
    theta_ep = pseudoeq_potential_T(T, P, w, P_0)
    return theta_ep
def equiv_potential_T(T, P, P_0 = 1000.):
    '''Takes temperature (C), pressure (hPa) and reference pressure (hPa) and returns equivalent potential temperature (K).  The calculation for mixing ratio (kg/kg) is included in this definition and need not be an input.'''
    e = 0.
    e_s = sat_vapor_pressure(T)
    w_s = sat_mixing_ratio(P, T)
    w_t = w_s 
    c_wd = c_p_dry + w_t*c_L
    theta_e = (T+C_to_K)*(P_0/(P-e))**(R_d/c_wd)*np.exp((L_v*w_t)/(c_wd*(T+C_to_K)))
    return theta_e
