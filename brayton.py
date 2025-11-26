import numpy as np
R_gas = 287.0  # J/kgK

def Brayton_IG(PR:float, cyc:dict) -> tuple:
    """define ideal gas Brayton cycle performance"""
    # return specific work, efficiency, entropy and temperature at each state point
    
    TR    = cyc['TR'].value       # temperature ratio (tmax/tmin)
    gam_c = cyc['gam_c'].value    # specific heat capacity ratio for compressor 
    gam_t = cyc['gam_t'].value    # specific heat capacity ratio for turbine 
    etap_c = cyc['etap_c'].value  # polytropic eff for compressor
    etap_t = cyc['etap_t'].value  # polytropic eff for turbine
    PR_cc  = cyc['PR_cc'].value   # combustion chamber pressure ratio

    
    # Inputs, gam_c = gamma_air, gam_t = gamma_flg
    T02 = 293.15
    cp_c = R_gas*gam_c/(gam_c-1)
    cp_t = R_gas*gam_t/(gam_t-1)
    
    # Compressor stage 
    TR_c = PR**((gam_c-1)/(gam_c*etap_c)) # Compressor temperature ratio
    T03 = T02*TR_c

    # Combustion 
    T04 = T02*TR

    # Pressure loss in combustion chamber 
    PR_t = PR*PR_cc

    # Turbine stage
    TR_t = (1/PR_t)**(((gam_t-1)/gam_t)*etap_t)
    T05 = TR_t*T04 # Essentially T01?
    
    # Cycle
    w_t = cp_t*(T04-T05)
    w_c = cp_c*(T03-T02)
    h_in = (cp_c+cp_t)/2*(T04-T03)
    eta = (w_t-w_c)/h_in
    specw = (w_t-w_c)/(cp_c*T02) # non-dimensional specific work of the cycle
    
    # Entropy and temperature
    ds_3 = cp_c*np.log(T03/T02) - R_gas*np.log(PR) # compression
    ds_4 = 0.5*(cp_c+cp_t)*np.log(T04/T03) - R_gas*np.log(PR_cc) # heat injection
    ds_5 = cp_t*np.log(T05/T04) - R_gas*np.log(1/PR_t) # expansion
    
    s_2 = 0
    s_3 = ds_3
    s_4 = s_3 + ds_4
    s_5 = s_4 + ds_5
    
    # Concatenate tuples
    entr = [s_2, s_3, s_4, s_5, s_2] # s_2 = 0, reference!
    temp = [T02, T03, T04, T05, T02]
    
    return specw, eta, entr, temp