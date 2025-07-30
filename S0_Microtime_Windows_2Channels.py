import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scripts.Read_PTU import read_data
from datetime import datetime
import json

# PARAMETER
#########################################################################################
DATA_FOLDER = "DDX3X"

DONOR_CHANNEL = 1 # detection channel of the donor signal
ACCEPTOR_CHANNEL = 2 # detection channel of the acceptor signal

NUM_CHANNELS = 0 # number of microtime channels - gets calculated from data files if zero
#########################################################################################

mpl.use('TkAgg') # uses external plotting

# search for HT3-files
contPath = [f for f in os.listdir(DATA_FOLDER) if f.endswith('.ptu')]

# read first file
data, unit, globRes, binRes = read_data(DATA_FOLDER + '/' + contPath[0], 1)

dt = np.round(binRes*1e12).astype(int) # (ps) temporal resolution of microtime channels

if NUM_CHANNELS == 0:

    NUM_CHANNELS = np.round(globRes/binRes).astype(int)

# separate channels
microD = data[data[:, 0] == DONOR_CHANNEL, 1]
microA = data[data[:, 0] == ACCEPTOR_CHANNEL, 1]

edges = np.arange(1, NUM_CHANNELS)

# calculate microtime histograms
hD, _ = np.histogram(microD, bins=np.append(edges, edges[-1] + 1))
hA, _ = np.histogram(microA, bins=np.append(edges, edges[-1] + 1))


# Histogram plotting
f1, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4))

# Donor - channel
ax1.semilogy(edges, hD, color=(0, 0.8, 0))
ax1.set_xlim(0, NUM_CHANNELS)
ax1.set_ylim(0.1, np.max(hD))
ax1.set_ylabel('Detector ' + str(DONOR_CHANNEL))
ylim1 = ax1.get_ylim()
ax1.set_ylim(1, ylim1[1])

# Acceptor - channel
ax2.semilogy(edges, hA, color=(0.8, 0, 0))
ax2.set_xlim(0, NUM_CHANNELS)
ax2.set_ylabel('Detector ' + str(ACCEPTOR_CHANNEL))
ax2.set_xlabel('Channel')
ylim2 = ax2.get_ylim()
ax2.set_ylim(1, ylim2[1])

plt.tight_layout()

# Range selection
print("Select FRET region (2 clicks)")
xFRET = plt.ginput(2)
xFRET = [x[0] for x in xFRET]

print(xFRET)

# Redraw with fill
ax1.fill_between([xFRET[0], xFRET[1]], ylim1[0], ylim1[1], color=(0.75, 0.75, 0.75))
ax1.semilogy(edges, hD, color=(0, 0.8, 0))
ax1.set_xlim(0, NUM_CHANNELS)
ax1.set_ylabel('Detector ' + str(DONOR_CHANNEL))

ax2.fill_between([xFRET[0], xFRET[1]], ylim2[0], ylim2[1], color=(0.75, 0.75, 0.75))
ax2.semilogy(edges, hA, color=(0.8, 0, 0))
ax2.set_xlim(0, NUM_CHANNELS)
ax2.set_ylabel('Detector ' + str(ACCEPTOR_CHANNEL))
ax2.set_xlabel('Channel')

print("Select Acceptor check region (2 clicks)")
xACC = plt.ginput(2)
xACC = [x[0] for x in xACC]

print(xACC)

# Fill both regions
ax2.fill_between([xACC[0], xACC[1]], ylim2[0], ylim2[1], color=(0.75, 0.75, 0.75))
ax2.fill_between([xFRET[0], xFRET[1]], ylim2[0], ylim2[1], color=(0.75, 0.75, 0.75))
ax2.semilogy(edges, hA, color=(0.8, 0, 0))
ax2.set_xlim(0, NUM_CHANNELS)
ax2.set_ylabel('Detector ' + str(ACCEPTOR_CHANNEL))
ax2.set_xlabel('Channel')

# Save settings
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
figpath = os.path.join("settings", f"Selected_time_windows_{timestamp}.png")
f1.savefig(figpath)

dt = np.round(binRes*1e12).astype(int)

settings_arr = np.round([*xFRET, *xACC]).astype(int)
settings_dict = {'FRET': settings_arr[:2].tolist(), 'Acceptor': settings_arr[2:].tolist(), 'Channels': NUM_CHANNELS.tolist(), 'dt': dt.tolist(), 'Donor_channel': DONOR_CHANNEL, 'Acceptor_channel': ACCEPTOR_CHANNEL}

settings_path = os.path.join("settings", f"Settings_{timestamp}.json")

with open(settings_path, 'w') as f:
    json.dump(settings_dict, f, indent=4)

plt.show()

