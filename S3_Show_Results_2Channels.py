import os
import tkinter as tk
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from tkinter import filedialog, Tk, messagebox
import json
from scripts.Scatter2Density import Scatter2Density
import pandas as pd

# PARAMETER
#########################################################################
RESULTS_PATH = "results"
RESULTS_FILE = "Results_20250926_hpT5_100mM_NaCl_PTU"

PATH_OUT='plots_and_histograms'

boolPNG = 1 # export plots and histograms as PNG
boolSVG = 1 # export plots and histograms as SVG
boolEXCEL = 1 # export results as Excel sheet

# correction factors for stoichiometry and FRET efficiency
ALPHA = 0.0416
BETA = 0.0236
GAMMA = 0.6273

# donor only lifetime
TAU_D0 = 2.5 # (ns)

# common burst filter to filter for double labeled molecules
NUM_PH = np.array([0, 500]) # minimal and maximal number of total photons
BRD_S = np.array([0.2, 0.8]) # lower and upper threshold of stoichiometry filter
BRD_ALEX2CDE = np.array([-1, 12]) # lower and upper threshold of ALEX-2CDE filter

# special burst filter - usually kept open!!!
BRD_E = np.array([-0.1, 1.1]) # lower and upper threshold of stoichiometry filter
BRD_FRET2CDE = np.array([-1000, 1000]) # lower and upper threshold of FRET-2CDE filter
RATIO_NGNR = np.array([-10, 10]) # lower and upper threshold of NG/NR filter
BRD_TAU_D = np.array([-100, 100]) # lower and upper threshold of donor lifetime
BRD_TAU_A = np.array([-100, 100]) # lower and upper threshold of acceptor lifetime

# settings for FRET histogram binning
BIN_SIZE = 0.03
OFFSET = 0.017
#########################################################################

# dividing a by b without b=0 warning
def div_array(a, b):
    return np.divide(a, b, out=np.full_like(a, np.inf, dtype=float), where=b != 0)

# makes it possible to save np arrays in json format
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

mpl.use('TkAgg') # uses external plotting

# Load data
BACKGROUND_FILE = os.path.join(RESULTS_PATH, f"BG_{RESULTS_FILE[8:]}.json")
SETTINGS_FILE = os.path.join(RESULTS_PATH, f"Settings_{RESULTS_FILE[8:]}.json")

# open data file
with open(os.path.join(RESULTS_PATH, f"{RESULTS_FILE}.json"), 'r') as f:
    results = json.load(f)

# open corresponding background file
with open(BACKGROUND_FILE, 'r') as f:
    background = json.load(f)

# open settings file
with open(SETTINGS_FILE, 'r') as f:
        settings = json.load(f)

# load parameters from result file
BN = np.array(results['BN']) # burst number
PosT = np.array(results['PosT']) # (s)
BIN_T = np.array(results['BIN_T']) # (ms)
ID = np.array(results['ID']) # (kHz)
IA = np.array(results['IA']) # (kHz)
IA0 = np.array(results['IA0']) # (kHz)
TauD = np.array(results['TauD']) # (ns)
TauA0 = np.array(results['TauA0']) # (ns)
FRET2CDE = np.array(results['FRET2CDE'])
ALEX2CDE = np.array(results['ALEX2CDE'])
DTGR_TR0 = np.array(results['DTGR_TR0']) # (ms)

if not BN.any():

    print('Results file is empty!')

else:

    # load background data
    BD = background['BD_mean'] # (kHz)
    BA = background['BA_mean'] # (kHz)
    BA0 = background['BA0_mean'] # (kHz)

    Nph = ID + IA + IA0 # (kHz) total number of photons

    # background corrected fluorescence intensity
    FD = ID - BIN_T * BD # (kHz)
    FA = IA - BIN_T * BA # (kHz)
    FA0 = IA0 - BIN_T * BA0 # (kHz)

    ratioDA = np.log10(div_array(ID, IA) + 1e-10) # intensity ratio

    # corrected FRET efficiency and stoichiometry
    E = div_array((FA - ALPHA * FA0 - BETA * FD), (FA - ALPHA * FA0 - BETA * FD + GAMMA * FD))
    S = div_array((FA - ALPHA * FA0 - BETA * FD + GAMMA * FD), (FA - ALPHA * FA0 - BETA * FD + GAMMA * FD + FA0))

    # relative fluorescence lifetime of the donor
    Rel_TauD = TauD / TAU_D0 # relative fluorescence lifetime of the donor

    idx = (BRD_E[0] <= E) & (E <= BRD_E[1]) \
          & (NUM_PH[0] <= Nph) & (Nph <= NUM_PH[1]) \
          & (BRD_S[0] <= S) & (S <= BRD_S[1]) \
          & (BRD_ALEX2CDE[0] <= ALEX2CDE) & (ALEX2CDE <= BRD_ALEX2CDE[1]) \
          & (BRD_FRET2CDE[0] <= FRET2CDE) & (FRET2CDE <= BRD_FRET2CDE[1]) \
          & (RATIO_NGNR[0] <= ratioDA) & (ratioDA <= RATIO_NGNR[1]) \
          & (BRD_TAU_D[0] <= TauD) & (TauD <= BRD_TAU_D[1]) \
          & (BRD_TAU_A[0] <= TauA0) & (TauA0 <= BRD_TAU_A[1])

    # data filtering
    sel_E = E[idx]
    sel_S = S[idx]
    sel_TauD = TauD[idx] # (ns)
    sel_Rel_TauD = Rel_TauD[idx]
    sel_TauA0 = TauA0[idx] # (ns)
    sel_ALEX2CDE = ALEX2CDE[idx]
    sel_FRET2CDE = FRET2CDE[idx]
    sel_PosT = PosT[idx] # (s)

    # edges for FRET efficiency histogram
    edgesE = np.arange(-0.1 + OFFSET, 1.1 + OFFSET, BIN_SIZE)
    hist_E, _ = np.histogram(sel_E, bins=np.append(edgesE, edgesE[-1] + BIN_SIZE))

    # plotting of results
    f1, axs = plt.subplots(2, 3, figsize=(16, 9), gridspec_kw={'wspace': 0.3, 'hspace': 0.3})
    ax1, ax2, ax3, ax4, ax5, ax6 = axs.flatten()
    f1.suptitle(RESULTS_FILE)

    # Create a custom colormap: white -> red
    cmap = LinearSegmentedColormap.from_list('white_red', ['white', 'red'])

    # FRET efficiency vs. Stoichiometry plot
    sel_E1, sel_S1, z_ES1 = Scatter2Density(sel_E, sel_S)

    ax1.scatter(E, S, c='black', s=2)
    ax1.scatter(sel_E1, sel_S1, c=z_ES1, s=7, cmap='jet')
    ax1.set_xlim(-0.1, 1.1)
    ax1.set_ylim(-0.1, 1.1)
    ax1.set_xlabel('FRET efficiency, $E$')
    ax1.set_ylabel('Stoichiometry, $S$')

    # FRET efficiency vs. FRET-2CDE plot
    sel_E2, sel_ALEX2CDE2, z_EA2E2 = Scatter2Density(sel_E, sel_ALEX2CDE)

    ax2.scatter(E, ALEX2CDE, c='black', s=2)
    ax2.scatter(sel_E2, sel_ALEX2CDE2, c=z_EA2E2, s=7, cmap='jet')
    ax2.set_xlim(-0.1, 1.1)
    ax2.set_ylim(0, 110)
    ax2.set_xlabel('FRET efficiency, $E$')
    ax2.set_ylabel('ALEX-2CDE')

    # FRET efficiency vs. FRET-2CDE plot
    sel_E3, sel_FRET2CDE3, z_EF2E3 = Scatter2Density(sel_E, sel_FRET2CDE)

    ax3.scatter(sel_E3, sel_FRET2CDE3, c=z_EF2E3, s=7, cmap='jet')
    ax3.set_xlim(-0.1, 1.1)
    ax3.set_ylim(0, 110)
    ax3.set_xlabel('FRET efficiency, $E$')
    ax3.set_ylabel('FRET-2CDE')

    # FRET efficiency histogram
    ax4.hist(sel_E, bins=np.append(edgesE, edgesE[-1] + BIN_SIZE), color='skyblue', edgecolor='black')
    ax4.set_xlim(-0.1, 1.1)
    ax4.set_xlabel('FRET efficiency, $E$')
    ax4.set_ylabel('Number of molecules')

    # FRET efficiency vs. relative fluorescence lifetime of the donor
    idx_number=~np.isnan(sel_Rel_TauD)
    sel_En = sel_E[idx_number]
    sel_Rel_TauDn = sel_Rel_TauD[idx_number] # NaN filtered molecules

    sel_En5, sel_Rel_TauDn5, z_ETD5 = Scatter2Density(sel_En, sel_Rel_TauDn)

    ax5.scatter(sel_En5, sel_Rel_TauDn5, c=z_ETD5, s=7, cmap='jet')
    ax5.plot(np.array([0, 1]), np.array([1, 0]), color='black')
    ax5.set_xlim(-0.1, 1.1)
    ax5.set_ylim(-0.1, 1.1)
    ax5.set_xlabel('FRET efficiency, $E$')
    ax5.set_ylabel(r'$\tau_D / \tau_{D0}$')

    # Fluorescence lifetime plots
    edgesT = np.arange(0, 10, 0.1)
    hist_TD_only, _ = np.histogram(TauD[S > 0.9], bins=np.append(edgesT, edgesT[-1] + 0.1))
    hist_TA_only, _ = np.histogram(TauA0[S < 0.2], bins=np.append(edgesT, edgesT[-1] + 0.1))
    hist_TD, _ = np.histogram(sel_TauD, bins=np.append(edgesT, edgesT[-1] + 0.1))

    ax6.step(edgesT, hist_TD_only / np.max(hist_TD_only), where='post', color=(0, 0.75, 0))
    ax6.step(edgesT, hist_TD / np.max(hist_TD), where='post', color=(0.75, 0.75, 0))
    ax6.step(edgesT, hist_TA_only / np.max(hist_TA_only), where='post', color=(0.75, 0, 0))
    ax6.set_xlim(0, 10)
    ax6.set_xlabel(r'$\tau$ (ns)')
    ax6.set_ylabel('number of events')
    ax6.legend(['Donor only', 'Donor - FRET', 'Acceptor only'])

    plt.show()

    root = tk.Tk()
    root.withdraw()  # Hide the main window

    # ask for saving data
    result = messagebox.askyesno("Confirm", "Do you want to save the results?")

    if result:

        # save parameters in json file
        if not os.path.exists(PATH_OUT):
            os.makedirs(PATH_OUT)

        PATH_OUT_FOLDER = os.path.join(PATH_OUT, RESULTS_FILE[8:])

        if not os.path.exists(PATH_OUT_FOLDER):
            os.makedirs(PATH_OUT_FOLDER)

        # extract and save plots
        bbox1 = ax1.get_tightbbox(f1.canvas.get_renderer()).expanded(1.025, 1.025)
        bbox2 = ax2.get_tightbbox(f1.canvas.get_renderer()).expanded(1.025, 1.025)
        bbox3 = ax3.get_tightbbox(f1.canvas.get_renderer()).expanded(1.025, 1.025)
        bbox4 = ax4.get_tightbbox(f1.canvas.get_renderer()).expanded(1.025, 1.025)
        bbox5 = ax5.get_tightbbox(f1.canvas.get_renderer()).expanded(1.025, 1.025)
        bbox6 = ax6.get_tightbbox(f1.canvas.get_renderer()).expanded(1.025, 1.025)

        # transform to figure coordinate system
        bbox1 = bbox1.transformed(f1.dpi_scale_trans.inverted())
        bbox2 = bbox2.transformed(f1.dpi_scale_trans.inverted())
        bbox3 = bbox3.transformed(f1.dpi_scale_trans.inverted())
        bbox4 = bbox4.transformed(f1.dpi_scale_trans.inverted())
        bbox5 = bbox5.transformed(f1.dpi_scale_trans.inverted())
        bbox6 = bbox6.transformed(f1.dpi_scale_trans.inverted())

        if boolPNG:

            f1.savefig(os.path.join(PATH_OUT_FOLDER, "S_vs_E_scatter_plot.png"), bbox_inches=bbox1, dpi=300, format='png')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "ALEX2CDE_vs_E_scatter_plot.png"), bbox_inches=bbox2, dpi=300, format='png')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "FRET2CDE_vs_E_scatter_plot.png"), bbox_inches=bbox3, dpi=300, format='png')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "E_histogram.png"), bbox_inches=bbox4, dpi=300, format='png')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "RelTauD_vs_E_scatter_plot.png"), bbox_inches=bbox5, dpi=300, format='png')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "Lifetime_histograms.png"), bbox_inches=bbox6, dpi=300, format='png')

        if boolSVG:

            f1.savefig(os.path.join(PATH_OUT_FOLDER, "S_vs_E_scatter_plot.svg"), bbox_inches=bbox1, format='svg')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "ALEX2CDE_vs_E_scatter_plot.svg"), bbox_inches=bbox2, format='svg')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "FRET2CDE_vs_E_scatter_plot.svg"), bbox_inches=bbox3, format='svg')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "E_histogram.svg"), bbox_inches=bbox4, format='svg')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "RelTauD_vs_E_scatter_plot.svg"), bbox_inches=bbox5,format='svg')
            f1.savefig(os.path.join(PATH_OUT_FOLDER, "Lifetime_histograms.svg"), bbox_inches=bbox6, format='svg')

        if boolEXCEL:

            df = pd.DataFrame({
                "time (s)": sel_PosT,
                "E": sel_E,
                "S": sel_S,
                "Tau_D (ns)": sel_TauD,
                "Tau_A0 (ns)": sel_TauA0,
                "ALEX2CDE": sel_ALEX2CDE,
                "FRET2CDE": sel_FRET2CDE
            })

            df.to_excel(os.path.join(PATH_OUT_FOLDER, "Results.xlsx"), index=False)

        if (boolPNG | boolSVG | boolEXCEL):

            # generate dictionary from data arrays
            corr_filter_dict = {

                'alpha': np.array([ALPHA]),
                'beta' : np.array([BETA]),
                'gamma' : np.array([GAMMA]),
                'TauD0' : np.array([TAU_D0]),
                'number_of_Photons' :NUM_PH,
                'borders_S' : BRD_S,
                'borders_ALEX2CDE' : BRD_ALEX2CDE,
                'borders_E': BRD_E,
                'borders_FRET2CDE': BRD_FRET2CDE,
                'borders_NGNR': RATIO_NGNR,
                'borders_TauD': BRD_TAU_D,
                'borders_TauA': BRD_TAU_A,
                'Bin_Size' : np.array([BIN_SIZE]),
                'Offset' : np.array([OFFSET])
            }

            corr_filter_path = os.path.join(PATH_OUT_FOLDER, "CorrFilter_settings.json")
            settings_path = os.path.join(PATH_OUT_FOLDER, "Analysis_settings.json")

            with open(corr_filter_path, 'w') as f:
                json.dump(corr_filter_dict, f, indent=4, cls=NumpyEncoder)

            with open(settings_path, 'w') as f:
                json.dump(settings, f, indent=4)
