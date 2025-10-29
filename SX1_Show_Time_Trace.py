import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scripts.Read_PTU import read_data
import json

# PARAMETER
#########################################################################
DATA_FOLDER = "20240725HFRuler"

SETTINGS_FILE = "Settings_20250930_215355"

BIN_T = 1  # (ms) bin time to calculate time trace
#########################################################################

mpl.use('TkAgg') # uses external plotting

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

data, unit, globRes, binRes = read_data(DATA_FOLDER + '/' + contPath[0])

# photon macrotimes
macroAll = data[(((BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1])) | ((BRD_ACC[0] < data[:, 1]) & (data[:, 1] < BRD_ACC[1]))), 2] * 1e-6 # (ms)

macroD = data[(BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1]) & (data[:, 0] == DONOR_CHANNEL), 2] * 1e-6 # (ms)
macroA = data[(BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1]) & (data[:, 0] == ACCEPTOR_CHANNEL), 2] * 1e-6 # (ms)
macroA0 = data[(BRD_ACC[0] < data[:, 1]) & (data[:, 1] < BRD_ACC[1]) & (data[:, 0] == ACCEPTOR_CHANNEL), 2] * 1e-6  # (ms)

# signal intensities
edges = np.arange(macroAll[0], macroAll[-1] + BIN_T, BIN_T) # (ms)

I_D, _ = np.histogram(macroD, bins=edges) # (kHz)
I_A, _ = np.histogram(macroA, bins=edges) # (kHz)
I_A0, _ = np.histogram(macroA0, bins=edges)  # (kHz)

I_All = I_D + I_A + I_A0

# start collection of background regions
fig = plt.figure(figsize=(10, 4))
gs = gridspec.GridSpec(3, 2)  # 1 row, 2 columns

ax1 = fig.add_subplot(gs[0, :])
ax2 = fig.add_subplot(gs[1, :])
ax3 = fig.add_subplot(gs[2, :])

ax1.plot(edges[1:], I_D, color=(0, 0.75, 0))
ax1.plot(edges[1:], -I_A, color=(0.75, 0, 0))
ax1.set_ylim(-1.1*np.max(I_A), 1.1*np.max(I_D))
ax1.set_xticklabels([])
ax1.set_ylabel('I_D & I_A (kHz)')

ax2.plot(edges[1:], I_A0, color=(0.75, 0, 0))
ax2.set_ylim(0, 1.1*np.max(I_A0))
ax2.set_xticklabels([])
ax2.set_ylabel('I_A0 (kHz)')

ax3.plot(edges[1:], I_All, color=(0, 0, 0))
ax3.set_ylim(0, 1.1*np.max(I_All))
ax3.set_ylabel('I_tot (kHz)')
ax3.set_xlabel('macrotime (ms)')

ax1.sharex(ax2)
ax2.sharex(ax1)
ax3.sharex(ax1)

plt.show()
