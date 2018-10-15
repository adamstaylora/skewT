import numpy as np


# constants

C_to_K = 273.15
c_P_dry = 1005.7
c_v_dry = 718
eps = 0.6220
k_dry = 0.2854


def sat_vapor_pressure(T):
    e_s = 6.112*np.exp((17.67*(T+C_to_K)/((T+C_to_K)+243.5)))
    return e_s
def sat_vapor_temperature(e_s):
    T = (243.5*np.log(e_s)-440.8)/(19.48-np.log(e_s))
    return T
def sat_mixing_ratio(P, T):
    w = eps*sat_vapor_pressure(T)/(P*100-sat_vapor_pressure(T))
    return w
def mixing_ratio_line(P, w_s):
    e_s = (w_s*P*100)/(w_s+eps)
    T = sat_vapor_temperature(e_s)
    return T
def RH(T, P, w):
    RH = 100*w/sat_mixing_ratio(P, T)
    return RH
def T_LCL(T, RH):
    T = 1/((1/(T+C_to_K-55))-(np.log(RH/100)/2840.))
    return T
def theta_dry(theta, P, P_0 = 1000.):
    T = (theta+C_to_K)*(P_0/P)**(-1*k_dry)
    return T
def pseudoeq_potential_T(T, P, w, P_0 = 1000.):
    theta_ep = (T+C_to_K)*(P_0/P)**(0.2854(1-0.28*w))
    return theta_ep
def theta_ep_field(T, P, P_0 = 1000.):
    theta_ep = pseudoeq_potential_T(T_LCL(T,RH(T,P,sat_mixing_ratio(P,T))), P, sat_mixing_ratio(P,T), P_0)
    return theta_ep

