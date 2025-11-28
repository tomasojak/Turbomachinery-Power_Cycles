import numpy as np
import matplotlib.pyplot as plt


from brayton import Brayton_IG
from IntercooledRecoup import RecoupIntercool_IG

defaultFontSize = 12


# fake a slider input
class property:
    value = None
    def __init__(self, value):
        self.value = value

class plotFunc:
    def __init__(self, func:callable, name:str):
        self.func = func
        self.name = name

# Because Python passes references, we need to reset the cycle dict each time we call this function
def resetCycle(cycle:dict, defaultCycle:dict) -> None:
    """Reset cycle dict to default values"""
    for key in defaultCycle:
        cycle[key].value = defaultCycle[key].value

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

def plotCycle():
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


def cycleComp():
    
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

def effiPerfPlot():
    
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

def cyclesEffies():
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

def prTrPerfPlot():
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

def etaWorkLoci():
     
    PRRange = np.linspace(1, 40, 10)
    TRRange = np.linspace(3, 6, 4)

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = ['tab:blue', 'tab:green', 'tab:orange', 'tab:red'] # abandnedoned

    for TR in TRRange:
        for PR in PRRange:
            resetCycle(cycle, defaultCycle)
            cycle['PR'].value = PR
            cycle['TR'].value = TR
            specw, eta, _, _ = RecoupIntercool_IG(cycle['PR'].value, cycle)
            ax.plot(specw, eta, 'o', color='black', markersize=4) # B&W version
            # ax.plot(specw, eta, 'o', color=colors[list(TRRange).index(TR)], markersize=4) # colored version
        # Ugly way to align
        if TR != TRRange[-1]:
            ax.text(specw + 0.15, eta, f'TR={TR:.1f}', fontsize=defaultFontSize)
        else:
            ax.text(specw - 0.1, eta - 0.05, f'TR={TR:.1f}', fontsize=defaultFontSize)

    ax.annotate('', xy=(1.05, 0.55), xytext=(0.65, 0.5),
                arrowprops=dict(arrowstyle='->', alpha=0.4, color='k'))
    text1 = ax.text(0.8, 0.495, 'Increasing TR', fontsize=defaultFontSize, backgroundcolor = 'white')

    ax.annotate('', xy=(1.3, 0.56), xytext=(1.1, 0.66),
                arrowprops=dict(arrowstyle='->', alpha=0.4, color='k'))
    text1 = ax.text(1.1, 0.65, 'Increasing PR', fontsize=defaultFontSize, backgroundcolor = 'white')
    ax.set_xlabel("Specific work (J/kg)")
    ax.set_ylabel("Efficiency (-)")

    ax.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.5)


##### MAIN #####

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
    S = np.array(entr)
    T = np.array(temp)

    ##### PLOTTING #####

    plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": defaultFontSize,
    "axes.labelsize": defaultFontSize,
    "axes.titlesize": defaultFontSize,
    "legend.fontsize": defaultFontSize,
    "xtick.labelsize": defaultFontSize,
    "ytick.labelsize": defaultFontSize,
    })

    plots = [
        plotFunc(plotCycle, "cycle_T-s_diagram"),
        plotFunc(cycleComp, "the_two_cycles_compared"),
        plotFunc(effiPerfPlot, "efficiencies_vs_performance"),
        plotFunc(cyclesEffies, "cycles_and_efficiencies"),
        plotFunc(prTrPerfPlot, "PR_and_TR_vs_performance"),
        plotFunc(etaWorkLoci, "eta_work_loci")
    ]

    savePlots = True
    saveDir = "./figures/"

    for plot in plots:
        plot.func()
        fig = plt.gcf()
        if savePlots:
            savepath = f"{saveDir}{plot.name}.pdf"
            plt.savefig(savepath, format='pdf', bbox_inches='tight')
            plt.close(fig)
            print(f"Saved plot {plot.name} to {savepath}")

    plt.show()