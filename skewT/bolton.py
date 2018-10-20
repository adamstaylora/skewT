import numpy as np


# constants

C_to_K = 273.15
c_P_dry = 1005.7
c_v_dry = 718
eps = 0.6220
k_dry = 0.2854


def sat_vapor_pressure(T):
    e_s = 6.112*np.exp((17.67*T)/(T+243.5)) #celsius input, hPa output
    return e_s
def sat_vapor_temperature(e_s):
    T = (243.5*np.log(e_s)-440.8)/(19.48-np.log(e_s)) #hPa input, celsius output
    return T
def sat_mixing_ratio(P, T):
    w = eps*sat_vapor_pressure(T)/(P-sat_vapor_pressure(T)) #hPa and celsius input, hPa from sat_vapor_pressure, kg/kg output
    return w
def mixing_ratio_line(P, w_s):
    e_s = (w_s*P)/(w_s+eps) #hPa and kg/kg input, hPa output
    T = sat_vapor_temperature(e_s) #celsius output
    return T
def RH(T, P, w):
    RH = 100*w/sat_mixing_ratio(P, T) #celsius and hPa input, kg/kg input, RH in %
    return RH
def T_LCL(T, RH):
    A = 1/(T+ C_to_K-55) # celsius input, RH % input
    B = np.log(RH/100.0)/2840.0
    T = 1/(A-B) #kelvin output
    return T
def theta_dry(theta, P, P_0 = 1000.):
    T = theta*(P_0/P)**(-1*k_dry) #kelvin input, kelvin output
    return T
def pseudoeq_potential_T(T, P, w, P_0 = 1000.0):
    rh = RH(T, P, w) #celsius, hPa, kg/kg input and kelvin output
    T_lcl = T_LCL(T, rh)
    A = ((3.376/T_lcl)-0.00254)*w*1000*(1+0.81*(w))
    theta_ep = np.exp(A)*(T+C_to_K)*(P_0/P)**(0.2854*(1-0.28*(w)))
    return theta_ep
def theta_ep_field(T, P, P_0 = 1000.0):
    w = sat_mixing_ratio(P,T)
    theta_ep = pseudoeq_potential_T(T, P, w, P_0)
    return theta_ep
