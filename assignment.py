import numpy as np

from brayton import Brayton_IG

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

    # TODO: Experiment with intercooler and recuperator efficiency values
    
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

# fake a slider input
class property:
    value = None
    def __init__(self, value):
        self.value = value

def resetCycle(cycle:dict, defaultCycle:dict) -> None:
    """Reset cycle dict to default values"""
    for key in defaultCycle:
        cycle[key].value = defaultCycle[key].value

# Because Python passes references, we need to reset the cycle dict each time we call this function
def valueRangeCycle(cycle:dict, cycleFunction:callable, targetX:str, Xvalues:list) -> tuple:
    """Evaluate cycle performance over a range of input values for a given property"""
    specw = np.zeros(len(Xvalues))
    eta = np.zeros(len(Xvalues))
    orig_val = cycle[targetX].value
    try:
        for iVal in range(len(Xvalues)):
            cycle[targetX].value = Xvalues[iVal]
            specw[iVal], eta[iVal], _, _ = cycleFunction(cycle['PR'].value, cycle)
    finally:
        cycle[targetX].value = orig_val
        resetCycle(cycle, defaultCycle)
    return specw, eta

if __name__ == "__main__":
    # Example usage
    defaultCycle = {
        'PR': property(32.0),
        'TR': property(6.0),
        'gam_c': property(1.4),
        'gam_t': property(1.4),
        'etap_c': property(1),
        'etap_t': property(1),
        'PR_cc': property(1)
    }
    cycle = {}
    for key in defaultCycle:
        cycle[key] = property(defaultCycle[key].value)
    
    specw, eta, entr, temp = RecoupIntercool_IG(cycle['PR'].value, cycle)
    print(f"Specific Work: {specw}, Efficiency: {eta}")
    print(f"Entropy states: {entr}")
    print(f"Temperature states: {temp}")

    # Plot T-s diagram
    import matplotlib.pyplot as plt
    S = np.array(entr)
    T = np.array(temp)

    # plt.figure(figsize=(6,5))
    # plt.semilogy(S, T, '-o', color='tab:tab:blue')  # log-scale for temperature (y-axis)
    # for i,(si,ti) in enumerate(zip(S,T)):
    #     plt.text(si, ti, f'  {i+2}', va='bottom', fontsize=9)  # label points 2..8
    # plt.xlabel('Entropy (J/kg/K)')
    # plt.ylabel('Temperature (K) [log scale]')
    # plt.title('T-s diagram: Recoup + Intercool Brayton (ideal gas)')
    # plt.grid(True, which='both', ls='--')
    # plt.gca().invert_xaxis() if False else None  # keep default orientation
    # plt.tight_layout()



    # fig, ax1 = plt.subplots()
    # ax2 = ax1.twinx()
    # ax1.plot(Xvalues, specw, 'g-')
    # ax2.plot(Xvalues, eta, 'b-')
    # ax1.set_ylabel('specw', color  ='g')
    # ax2.set_ylabel('eta', color = 'b')
    # ax1.set_xlabel(targetX)

    # plt.grid(True, which='both', ls='--')

    ##### PLOTTING #####

    fontSize = 10

    plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": fontSize,
    "axes.labelsize": fontSize,
    "axes.titlesize": fontSize,
    "legend.fontsize": fontSize,
    "xtick.labelsize": fontSize,
    "ytick.labelsize": fontSize,
    })

    # region ########### Plot cycle T-s ###########

    _, _, entrRecoup, tempRecoup = RecoupIntercool_IG(cycle['PR'].value, cycle)
    _, _, entrBrayton, tempBrayton = Brayton_IG(cycle['PR'].value, cycle)

    S_recoup = np.array(entrRecoup)
    T_recoup = np.array(tempRecoup)

    S_brayton = np.array(entrBrayton)
    T_brayton = np.array(tempBrayton)

    fig, ax = plt.subplots(figsize=(6,5))
    ax.semilogy(S_brayton, T_brayton, '-o', color='tab:blue', label='Brayton')  # log-scale for temperature (y-axis)
    ax.semilogy(S_recoup, T_recoup, '-o', color='tab:orange', label='Recoup + Intercool')  # log-scale for temperature (y-axis) 
    # for i,(si,ti) in enumerate(zip(S_brayton[:-1], T_brayton[:-1])):  # skip last point
    #     ax.text(si, ti, f'  {i+2}', va='bottom', fontsize=9)  # label points 2..(n-1)
    for i,(si,ti) in enumerate(zip(S_recoup[:-1], T_recoup[:-1])):  # skip last point
        ax.text(si, ti, f'  {i+2}', va='bottom', fontsize=9)  # label points 2..(m-1)

    ax.set_xlabel('Entropy (J/kg/K)')
    ax.set_ylabel('Temperature (K) [log scale]')
    
    ax.set_title('T-s diagram: Brayton vs Recoup + Intercool (ideal gas)')
    ax.grid(True, which='both', ls='--')
    ax.legend()

    plt.tight_layout()
    
    # endregion

    # region ########### Plot 1: the two cycles compared ###########

    Xvalues = np.linspace(1, 40, 50)
    
    specwBraytonIG, etaBraytonIG = valueRangeCycle(cycle, Brayton_IG, 'PR', Xvalues)
    specwRecoupICIG, etaRecoupICIG = valueRangeCycle(cycle, RecoupIntercool_IG, 'PR', Xvalues)

    fig, ax1 = plt.subplots(figsize=(8, 5))

    # -----------------------
    #   Primary y-axis: Spec. work
    # -----------------------

    ax1.plot(Xvalues, specwBraytonIG, linewidth=1,
            label="Brayton", color='tab:blue')
    ax1.plot(Xvalues, specwRecoupICIG, linewidth=1,
            label="Recoup. Intercooled", color='tab:blue', linestyle="--")

    ax1.set_xlabel("Pressure Ratio (-)")
    ax1.set_ylabel("Specific work (J/kg)", color='tab:blue')
    ax1.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.5)

    # -----------------------
    #   Vertical reference line
    # -----------------------
    PR_ref = 32
    ax1.axvline(PR_ref, color="black", linestyle="--", linewidth=1)
    ax1.text(PR_ref + 0.5, ax1.get_ylim()[1]*0.6, "PR = 32",
            ha="left", va="top")

    # -----------------------
    #   Secondary y-axis: Efficiency
    # -----------------------
    ax2 = ax1.twinx()

    ax2.plot(Xvalues, etaBraytonIG, linewidth=1,
            label="Brayton", color='tab:orange')
    ax2.plot(Xvalues, etaRecoupICIG, linestyle="--", linewidth=1,
            label="Recoup. Intercooled", color='tab:orange')

    ax2.set_ylabel("Efficiency (-)", color='tab:orange')

    # -----------------------
    #   Combined legend
    # -----------------------
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
            loc="lower right", frameon=True, framealpha=1)
    
    resetCycle(cycle, defaultCycle)
    # endregion

    # region ########### Plot 2: efficiencies vs performance ###########

    Xvalues = np.linspace(0.4, 1.0, 50)

    specwRecoupICIG_c, etaRecoupICIG_c = valueRangeCycle(cycle, RecoupIntercool_IG, 'etap_c', Xvalues)
    specwRecoupICIG_t, etaRecoupICIG_t = valueRangeCycle(cycle, RecoupIntercool_IG, 'etap_t', Xvalues)

    fig, ax1 = plt.subplots(figsize=(8, 5))

    # -----------------------
    #   Primary y-axis: Spec. work
    # -----------------------

    ax1.plot(Xvalues, specwRecoupICIG_c, linewidth=1,
            label="$w_s \sim \eta_{c}$", color='orange')
    ax1.plot(Xvalues, specwRecoupICIG_t, linewidth=1, linestyle="--",
            label="$w_s \sim \eta_{t}$", color='orange')
    
    ax1.set_xlabel("(Compressor / Turbine) Efficiency (-)")
    ax1.set_ylabel("Specific work (J/kg s)", color='orange')
    ax1.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.5)

    # -----------------------
    #   Secondary y-axis: Efficiency
    # -----------------------

    ax2 = ax1.twinx()
    ax2.plot(Xvalues, etaRecoupICIG_c, linestyle="--", linewidth=1,
            label="$\eta \sim \eta_{c}$", color='tab:blue')
    ax2.plot(Xvalues, etaRecoupICIG_t, linewidth=1,
            label="$\eta \sim \eta_{t}$", color='tab:blue')
    
    ax2.set_ylabel("Efficiency (-)", color='tab:blue')

    # -----------------------
    #   Combined legend
    # -----------------------
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
            loc="lower right", frameon=True, framealpha=1)

    # endregion

    # region ########### Cycles and efficiencies ###########

    Xvalues = np.linspace(0.65, 1.0, 50)

    _, etaRecoupICIG_c = valueRangeCycle(cycle, RecoupIntercool_IG, 'etap_c', Xvalues)
    _, etaBraytonIG_c = valueRangeCycle(cycle, Brayton_IG, 'etap_c', Xvalues)

    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(Xvalues, etaRecoupICIG_c, linewidth=1,
            label="recoup + intercool", color='black')
    ax1.plot(Xvalues, etaBraytonIG_c, linewidth=1,
            label="Brayton", color='black', linestyle='--')
    
    ax1.set_xlabel("Compressor Efficiency (-)")
    ax1.set_ylabel("Efficiency (-)")

    ax1.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.5)

    ax1.legend(loc="lower right", frameon=True, framealpha=1)

    resetCycle(cycle, defaultCycle)
    # endregion

    # region ########### PR and TR vs performance ###########

    # PR: Pressure ratio sweep
    PRSweep = np.linspace(1, 40, 100)
    specwRecoupPR, etaRecoupPR = valueRangeCycle(cycle, RecoupIntercool_IG, 'PR', PRSweep)

    # TR: Temperature ratio sweep
    TRSweep = np.linspace(3, 5, 100)
    specwRecoupTR, etaRecoupTR = valueRangeCycle(cycle, RecoupIntercool_IG, 'TR', TRSweep)

    # Create two vertically stacked subplots
    fig, (ax_pr, ax_tr) = plt.subplots(2, 1, figsize=(8, 10), constrained_layout=True)

    # -----------------------
    # Top: PR vs specific work (left) and efficiency (right)
    # -----------------------
    ax1 = ax_pr
    ax1.plot(PRSweep, specwRecoupPR, label="spec. work", color='black', linewidth=1)
    ax1.set_xlabel("Pressure Ratio (-)")
    ax1.set_ylabel("Specific work (J/kg)")
    ax1.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.5)

    ax1b = ax1.twinx()
    ax1b.plot(PRSweep, etaRecoupPR, linestyle="--", color='black', label="efficiency", linewidth=1)
    ax1b.set_ylabel("Efficiency (-)")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="lower right", frameon=True, framealpha=1)

    resetCycle(cycle, defaultCycle)

    # -----------------------
    # Bottom: TR vs specific work (left) and efficiency (right)
    # -----------------------
    ax2 = ax_tr
    ax2.plot(TRSweep, specwRecoupTR, label="spec. work", color='black', linewidth=1)
    ax2.set_xlabel("Temperature Ratio (Tmax/Tmin) (-)")
    ax2.set_ylabel("Specific work (J/kg)")
    ax2.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.5)

    ax2b = ax2.twinx()
    ax2b.plot(TRSweep, etaRecoupTR, linestyle="--", color='black', label="efficiency", linewidth=1)
    ax2b.set_ylabel("Efficiency (-)")

    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="lower right", frameon=True, framealpha=1)

    plt.show()

    # endregion