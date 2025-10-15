import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scripts.Read_PTU import read_data
import json

# PARAMETER
#########################################################################
DATA_FOLDER = "20250926_hpT5_100mM_NaCl_PTU"

SETTINGS_FILE = "Settings_20251001_214343"

BIN_T = 1  # (ms) bin time to calculate time trace
FRAC_T = 0.5  # (s) trace part to be analyzed
NUM_F = 5 # number of files ...
REG_F = 10 # ... and number of regions per file for background estimation
#########################################################################

mpl.use('TkAgg') # uses external plotting

def mean_bin_counts(x1, x2, I, edges):
    return np.mean(I[(edges[:-1] > x1 * 1000) & (edges[:-1] < x2 * 1000)])

# Load settings
with open(os.path.join("settings", f"{SETTINGS_FILE}.json"), 'r') as f:
    settings = json.load(f)

BRD_FRET = settings['FRET'] # microtime borders for FRET
BRD_ACC = settings['Acceptor'] # microtime borders for acceptor check
DONOR_CHANNEL = settings['Donor_channel'] # channel of the donor signal
ACCEPTOR_CHANNEL = settings['Acceptor_channel'] # channel of the acceptor signal

# Load measurement folder
FOLDER = os.path.basename(DATA_FOLDER)

contPath = [f for f in os.listdir(DATA_FOLDER) if f.endswith('.ptu')]

# start collection of background regions
f1, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

arrBD, arrBA, arrBA0, arrL = [], [], [], []

REG_F = REG_F + 1

for iterF in range(NUM_F):

    idx = iterF * len(contPath) // NUM_F

    data, unit, globRes, binRes = read_data(DATA_FOLDER + '/' + contPath[idx])

    # photon macrotimes
    macroAll = data[(((BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1])) | ((BRD_ACC[0] < data[:, 1]) & (data[:, 1] < BRD_ACC[1]))), 2] * 1e-6 # (ms)
    macroD = data[(BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1]) & (data[:, 0] == DONOR_CHANNEL), 2] * 1e-6 # (ms)
    macroA = data[(BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1]) & (data[:, 0] == ACCEPTOR_CHANNEL), 2] * 1e-6 # (ms)
    macroA0 = data[(BRD_ACC[0] < data[:, 1]) & (data[:, 1] < BRD_ACC[1]) & (data[:, 0] == ACCEPTOR_CHANNEL), 2] * 1e-6  # (ms)

    # signal intensities
    edges = np.arange(macroAll[0], macroAll[-1] + BIN_T, BIN_T) # (ms)

    I_All, _ = np.histogram(macroAll, bins=edges) # (kHz)

    I_D, _ = np.histogram(macroD, bins=edges) # (kHz)
    I_A, _ = np.histogram(macroA, bins=edges) # (kHz)
    I_A0, _ = np.histogram(macroA0, bins=edges)  # (kHz)

    measT = macroAll[-1] # total measurement time in milliseconds
    numT = int(np.floor(measT / (FRAC_T * 1000))) # number of time windows of current trace

    for iterR in range(1, min(numT, REG_F)):

        sel_edges = edges[((iterR * FRAC_T * 1000) < edges) & (edges < ((iterR + 1) * FRAC_T * 1000))] / 1000 # (s)
        sel_I_All = I_All[((iterR * FRAC_T * 1000) < edges[0:-1]) & (edges[0:-1] < ((iterR + 1) * FRAC_T * 1000))] # (kHz)

        # manual pick of background region
        if iterR > 1:
            line.remove()

        ax1.cla()
        line, = ax1.plot(sel_edges, sel_I_All, color=(0, 0, 0.75))
        ax1.set_xlim(sel_edges[0], sel_edges[-1])
        ax1.set_ylim(0, 100)
        ax1.set_xlabel('macrotime (s)')
        ax1.set_ylabel('total intenisty (kHz)')

        # Range selection
        xB = plt.ginput(2)
        xB = [x[0] for x in xB]

        arrBD.append(mean_bin_counts(xB[0], xB[1], I_D, edges))
        arrBA.append(mean_bin_counts(xB[0], xB[1], I_A, edges))
        arrBA0.append(mean_bin_counts(xB[0], xB[1], I_A0, edges))
        arrL.append(xB[1] - xB[0])

        print(f'{iterR}/{min(numT, REG_F) - 1} regions')

    print(f'{iterF + 1}/{NUM_F} files processed')

    edgesBG_D = np.arange(0, max(arrBD) + 4*max(arrBD)/10, max(arrBD)/10)  # (kHz)
    edgesBG_A = np.arange(0, max(arrBA) + 4*max(arrBA)/10, max(arrBA) / 10)  # (kHz)
    edgesBG_A0 = np.arange(0, max(arrBA0) + 4*max(arrBA0)/10, max(arrBA0) / 10)  # (kHz)

    hBD, _ = np.histogram(arrBD, bins=edgesBG_D)  # (kHz)
    hBA, _ = np.histogram(arrBA, bins=edgesBG_A)  # (kHz)
    hBA0, _ = np.histogram(arrBA0, bins=edgesBG_A0)  # (kHz)

    ax2.cla()
    ax2.step(edgesBG_D[1:], hBD, where='post', color=(0, 0.75, 0))
    ax2.step(edgesBG_A[1:], hBA, where='post', color=(0.75, 0, 0))
    ax2.step(edgesBG_A0[1:], hBA0, where='post', color=(0.75, 0, 0.75))
    ax2.set_xlim(0, max(max(arrBD), max(arrBA), max(arrBA0)))
    ax2.set_xlabel('background (kHz)')
    ax2.set_ylabel('number of events')

arrBD = np.array(arrBD)
arrBA = np.array(arrBA)
arrBA0 = np.array(arrBA0)

BD_mean = np.average(arrBD, weights=arrL)
BA_mean = np.average(arrBA, weights=arrL)
BA0_mean = np.average(arrBA0, weights=arrL)

BD_std = np.sqrt(np.average((arrBD - BD_mean) ** 2, weights=arrL))
BA_std = np.sqrt(np.average((arrBA - BA_mean) ** 2, weights=arrL))
BA0_std = np.sqrt(np.average((arrBA0 - BA0_mean) ** 2, weights=arrL))

background_dict = {
    'BD_mean': BD_mean,
    'BD_std': BD_std,
    'BA_mean': BA_mean,
    'BA_std': BA_std,
    'BA0_mean': BA0_mean,
    'BA0_std': BA0_std
}

if not os.path.exists("results"):
    os.makedirs("results")

background_path = os.path.join("results", f"BG_{FOLDER}.json")

with open(background_path, 'w') as f:
    json.dump(background_dict, f, indent=4)

print(background_dict)

plt.show()