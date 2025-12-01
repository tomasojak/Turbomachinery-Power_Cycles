import numpy as np
R_gas = 287.0  # J/kgK

def RecoupIntercool_IG(PR:float, cyc:dict) -> tuple:

    #########
    # Brayton cycle with intercooling + recuperation
    #########
    
    # Note: The implementation should not know this is a jupyter slider...
    # PR    = cyc['PR'].value       # pressure ratio
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

    #------------------------------

    # LPC
    #   Isentropic compression
    # We want both our compressor stages to have the same pressure ratio
    # Logical, but also specified in assignment text
    PR_c1 = PR**0.5
    PR_c2 = PR_c1
    
    TR_c1 = PR_c1**((gam_c-1)/(gam_c*etap_c)) # Compressor 1 temperature ratio
    T03 = T02*TR_c1
    
    # Intercooler stage - Isobaric cooling
        # intercooler efficiency
        # Let us define the intercooler efficiency as how close to T02 we cool
    etap_ic = 1 # Ideal intercooler, specified in assignment text

    T04 = T03 - etap_ic*(T03 - T02)

    # HPC
    #   Isentropic compression
    TR_c2 = PR_c2**((gam_c-1)/(gam_c*etap_c)) # Compressor 2 temperature ratio
    T05 = T04*TR_c2

    # Combustion stage
    T06 = T02*TR

    # Pressure loss in combustion chamber
    PR_t = PR*PR_cc

    # Turbine stage
    TR_t = (1/PR_t)**(((gam_t-1)/gam_t)*etap_t)
    T07 = TR_t*T06

    # Cycle
    w_t = cp_t*(T06-T07)
    w_c1 = cp_c*(T03-T02)
    w_c2 = cp_c*(T05-T04)
    
    if T07 > T05:
        hin = (cp_c+cp_t)/2*(T06-T07)
    else:
        hin = (cp_c+cp_t)/2*(T06-T05)

    eta = (w_t-(w_c1+w_c2))/hin
    specw = (w_t-(w_c1+w_c2))/(cp_c*T02) # non-dimensional specific work of the cycle    

    # Entropy and temperature
    ds_3 = cp_c*np.log(T03/T02) - R_gas*np.log(PR_c1)   # compression 1
    ds_4 = cp_c*np.log(T04/T03)                         # intercooling
    ds_5 = cp_c*np.log(T05/T04) - R_gas*np.log(PR_c2)   # compression 2
    ds_6 = cp_t*np.log(T06/T05)                         # combustion
    ds_7 = cp_t*np.log(T07/T06) + R_gas*np.log(PR_t)    # expansion
    
    s_2 = 0
    s_3 = ds_3
    s_4 = s_3 + ds_4
    s_5 = s_4 + ds_5
    s_6 = s_5 + ds_6
    s_7 = s_6 + ds_7
    
    # Concatenate tuples
    entr = [s_2, s_3, s_4, s_5, s_6, s_7, s_2] # s_2 = 0, reference!
    temp = [T02, T03, T04, T05, T06, T07, T02]
    
    return specw, eta, entr, temp