import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scripts.Read_PTU import read_data
import json

# PARAMETER
#########################################################################
DATA_FOLDER = "IRF_files"
DATA_FILE_DONOR = "IRF_ATTO532_saturated_KI_530nm"
DATA_FILE_ACCEPTOR = "IRF_ATTO655_saturated_KI_640nm"

SETTINGS_FILE = "Settings_20251029_180431"
#########################################################################

mpl.use('TkAgg') # uses external plotting

# makes it possible to save np arrays in json format
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

# Load settings
with open(os.path.join("settings", f"{SETTINGS_FILE}.json"), 'r') as f:
    settings = json.load(f)

BRD_FRET = settings['FRET'] # microtime borders for FRET
BRD_ACC = settings['Acceptor'] # microtime borders for acceptor check
NUM_CH = settings['Channels'] # number of non-zero microtime channels
DT_BIN = settings['dt'] # (ps) microtime resolution per channel
DONOR_CHANNEL = settings['Donor_channel'] # channel of the donor signal
ACCEPTOR_CHANNEL = settings['Acceptor_channel'] # channel of the acceptor signal

# loading of the donor and acceptor IRF
dataD, unitD, globResD, binResD = read_data(DATA_FOLDER + '/' + DATA_FILE_DONOR + '.ptu')
dataA, unitA, globResA, binResA = read_data(DATA_FOLDER + '/' + DATA_FILE_ACCEPTOR + '.ptu')

# microtime channels
microD = dataD[(BRD_FRET[0] < dataD[:, 1]) & (dataD[:, 1] < BRD_FRET[1]) & (dataD[:, 0] == DONOR_CHANNEL), 1]  # channel
microA = dataA[(BRD_ACC[0] < dataA[:, 1]) & (dataA[:, 1] < BRD_ACC[1]) & (dataA[:, 0] == ACCEPTOR_CHANNEL), 1]  # channel

# microtime bin edges
edgesD = np.arange(BRD_FRET[0], BRD_FRET[1])
edgesA = np.arange(BRD_ACC[0], BRD_ACC[1])

# calculation of microtime histograms
hD, _ = np.histogram(microD, bins=np.append(edgesD, edgesD[-1] + 1))
hA, _ = np.histogram(microA, bins=np.append(edgesA, edgesA[-1] + 1))

# background substraction
hD_bg = hD - np.mean(hD[-100:])
hA_bg = hA - np.mean(hA[-100:])

# calculation of the mean delay time
MEAN_IRF_DONOR = (np.sum(hD_bg * edgesD) / np.sum(hD_bg)) * DT_BIN * 1e-3  # (ns)
MEAN_IRF_ACCEPTOR = (np.sum(hA_bg * edgesA) / np.sum(hA_bg)) * DT_BIN * 1e-3  # (ns)

# Histogram plotting
f1, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4))

# Donor - channel
ax1.semilogy(edgesD * DT_BIN * 1e-3, hD_bg, color=(0, 0.8, 0))
ylim1 = ax1.get_ylim()
ax1.plot([MEAN_IRF_DONOR, MEAN_IRF_DONOR], [ylim1[0], ylim1[1]], color=(0, 0.8, 0), linestyle='--', label=f"<t> = {MEAN_IRF_DONOR:.4f} ns")
ax1.set_xlim(edgesD[1] * DT_BIN * 1e-3, edgesD[-1] * DT_BIN * 1e-3)
ax1.set_ylabel('Donor channel')
ax1.legend()

# Acceptor - channel
ax2.semilogy(edgesA * DT_BIN * 1e-3, hA_bg, color=(0.8, 0, 0))
ylim2 = ax2.get_ylim()
ax2.plot([MEAN_IRF_ACCEPTOR, MEAN_IRF_ACCEPTOR], [ylim2[0], ylim2[1]], color=(0.8, 0, 0), linestyle='--', label=f"<t> = {MEAN_IRF_ACCEPTOR:.4f} ns")
ax2.set_xlim(edgesA[1] * DT_BIN * 1e-3, edgesA[-1] * DT_BIN * 1e-3)
ax2.set_ylabel('Acceptor channel')
ax2.set_xlabel('Time (ns)')
ax2.legend()

# Save graph in settings
figpath = os.path.join("settings", f"IRF_mean_delay_time_{SETTINGS_FILE[9:]}.png")
f1.savefig(figpath)

plt.show()