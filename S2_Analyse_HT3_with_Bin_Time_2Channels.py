import os
import numpy as np
import matplotlib as mpl
from scripts.To_CDE_Functions import FRET_2CDE, ALEX_2CDE
from scripts.Read_PTU import read_data
from alive_progress import alive_bar
import time
import json

# PARAMETER
#########################################################################
DATA_FOLDER = "20240725HFRuler"

SETTINGS_FILE = "Settings_20250730_094959"

BIN_T = 1 # (ms) bin time
THRE_B = 30 # (kHz) lower intensity threshold for bin selection see ALGORITHM
ALGORITHM = 0 # 0 -> I_All >= THRE_B ... total intensity filter
              # 1 -> ((I_D + I_A) >= THRE_B) & (I_A0  >= THRE_B) ... demanding both active donor and active acceptor
              # 2 -> I_A0  >= THRE_B ... demanding active acceptor
              # 3 -> (I_D + I_A)  >= THRE_B ... demanding active donor
              # 4 -> I_D  >= THRE_B ... demanding donor counts

MEAN_IRF_DONOR = 2 # (ns) mean lifetime of the donor IRF
MEAN_IRF_ACCEPTOR = 27 # (ns) mean lifetime of the acceptor IRF
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

MID_CH = round(NUM_CH/2)

# Load measurement folder
FOLDER = os.path.basename(DATA_FOLDER)

contPath = [f for f in os.listdir(DATA_FOLDER) if f.endswith('.ptu')]

numF = len(contPath)

# arrays to collect the data
data_BN = np.empty((0, 1))
data_PosT = np.empty((0, 1))
data_BIN_T = np.empty((0, 1))
data_ID = np.empty((0, 1))
data_IA = np.empty((0, 1))
data_IA0 = np.empty((0, 1))
data_TauD = np.empty((0, 1))
data_TauA0 = np.empty((0, 1))
data_FRET2CDE = np.empty((0, 1))
data_ALEX2CDE = np.empty((0, 1))
data_DTGR_TR0 = np.empty((0, 1))

BN = 1

print(' ')
print('=====================================================')
print('Burst analysis running...')

with alive_bar(numF, force_tty=True) as bar:

    for iterF in range(numF):

        data, unit, globRes, binRes = read_data(DATA_FOLDER + '/' + contPath[iterF], 0)

        # photon arrival times -> macrotimes
        macroAll = data[(((BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1])) | ((BRD_ACC[0] < data[:, 1]) & (data[:, 1] < BRD_ACC[1]))), 2] * 1e-6 # (ms)
        macroD = data[(BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1]) & (data[:, 0] == DONOR_CHANNEL), 2] * 1e-6 # (ms)
        macroA = data[(BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1]) & (data[:, 0] == ACCEPTOR_CHANNEL), 2] * 1e-6  # (ms)
        macroA0 = data[(BRD_ACC[0] < data[:, 1]) & (data[:, 1] < BRD_ACC[1]) & (data[:, 0] == ACCEPTOR_CHANNEL), 2] * 1e-6  # (ms)

        macroDA = data[(BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1]) & ((data[:, 0] == DONOR_CHANNEL) | (data[:, 0] == ACCEPTOR_CHANNEL)), 2] * 1e-6  # (ms)

        # photon delay times -> microtimes
        microD = data[(BRD_FRET[0] < data[:, 1]) & (data[:, 1] < BRD_FRET[1]) & (data[:, 0] == DONOR_CHANNEL), 1] * DT_BIN * 1e-3  # (ns)
        microA0 = data[(BRD_ACC[0] < data[:, 1]) & (data[:, 1] < BRD_ACC[1]) & (data[:, 0] == ACCEPTOR_CHANNEL), 1] * DT_BIN * 1e-3  # (ns)

        lenT = (macroAll[-1] - macroAll[0]) / 1000 # (s)

        edges = np.arange(0, macroAll[-1], BIN_T) # (ms)

        # calculate histograms
        I_D, _ = np.histogram(macroD, bins=np.append(edges, edges[-1] + BIN_T))  # (kHz)
        I_A, _ = np.histogram(macroA, bins=np.append(edges, edges[-1] + BIN_T))  # (kHz)
        I_A0, _ = np.histogram(macroA0, bins=np.append(edges, edges[-1] + BIN_T))  # (kHz)

        I_All = I_D + I_A + I_A0 # (kHz) total intensity

        # algorithm for threshold filter to identify fluorescence bursts
        if ALGORITHM == 0:

            idx = I_All >= THRE_B

        elif ALGORITHM == 1:

            idx = ((I_D + I_A) >= THRE_B) & (I_A0  >= THRE_B)

        elif ALGORITHM == 2:

            idx = I_A0  >= THRE_B

        elif ALGORITHM == 3:

            idx = (I_D + I_A)  >= THRE_B

        elif ALGORITHM == 4:

            idx = I_D  >= THRE_B

        posB = np.where(idx == 1)[0] # burst positions
        numB = len(posB) # number of bursts

        # collect microtime information inside bursts -> arrays of parameters
        arrTauD = np.zeros(numB)
        arrTauA0 = np.zeros(numB)
        arrFRET2CDE = np.zeros(numB)
        arrALEX2CDE = np.zeros(numB)
        arrDTGR_TR0 = np.zeros(numB)

        arrPosT = np.zeros(numB)

        for iterA in range(numB):

            sub_macroAll = macroAll[(edges[posB[iterA]] <= macroAll) & ( macroAll < (edges[posB[iterA]] + BIN_T))] # (ms)
            sub_macroD = macroD[(edges[posB[iterA]] <= macroD) & (macroD < (edges[posB[iterA]] + BIN_T))] # (ms)
            sub_macroA = macroA[(edges[posB[iterA]] <= macroA) & ( macroA < (edges[posB[iterA]] + BIN_T))] # (ms)
            sub_macroA0 = macroA0[(edges[posB[iterA]] <= macroA0) & (macroA0 < (edges[posB[iterA]] + BIN_T))]  # (ms)

            sub_macroDA = macroDA[(edges[posB[iterA]] <= macroDA) & (macroDA < (edges[posB[iterA]] + BIN_T))] # (ms)

            sub_microD = microD[(edges[posB[iterA]] <= macroD) & (macroD < (edges[posB[iterA]] + BIN_T))]  # (ns)
            sub_microA0 = microA0[(edges[posB[iterA]] <= macroA0) & (macroA0 < (edges[posB[iterA]] + BIN_T))]  # (ns)

            # fluorescence lifetime calculation
            if len(sub_microD) == 0:

                arrTauD[iterA] = np.nan
            else:
                arrTauD[iterA] = np.mean(sub_microD) - MEAN_IRF_DONOR  # (ns)

            if len(sub_microA0) == 0:

                arrTauA0[iterA] = np.nan
            else:
                arrTauA0[iterA] = np.mean(sub_microA0) - MEAN_IRF_ACCEPTOR  # (ns)

            # Photon density indicators
            if (len(sub_macroDA) == 0) or (len(sub_macroD) == 0):

                arrFRET2CDE[iterA] = 0

            else:
                arrFRET2CDE[iterA] = FRET_2CDE(sub_macroDA, sub_macroD, 0.045) # kernel size is taken from the paper

            # Photon density indicators
            if (len(sub_macroA0) == 0) or (len(sub_macroDA) == 0):

                arrALEX2CDE[iterA] = 100

            else:
                arrALEX2CDE[iterA] = ALEX_2CDE(sub_macroA0, sub_macroDA, 0.075) # kernel size is taken from the paper

            if (len(sub_macroA0) == 0):

                arrDTGR_TR0[iterA] = -9.9

            elif (len(sub_macroDA) == 0):

                arrDTGR_TR0[iterA] = 9.9
            else:
                TR0_ = np.sum(sub_macroA0) / len(sub_macroA0)
                TGR_ = np.sum(sub_macroDA) / len(sub_macroDA)
                arrDTGR_TR0[iterA] = TGR_ - TR0_

            arrPosT[iterA] = np.mean(sub_macroAll) / 1000 + iterF * lenT # (s)

        # add data to output arrays
        data_BN = np.append(data_BN, np.arange(BN, BN + numB))
        data_PosT = np.append(data_PosT, arrPosT)
        data_BIN_T = np.append(data_BIN_T, BIN_T * np.ones(numB))
        data_ID = np.append(data_ID, I_D[idx])
        data_IA = np.append(data_IA, I_A[idx])
        data_IA0 = np.append(data_IA0, I_A0[idx])
        data_TauD = np.append(data_TauD, arrTauD)
        data_TauA0 = np.append(data_TauA0, arrTauA0)
        data_FRET2CDE = np.append(data_FRET2CDE, arrFRET2CDE)
        data_ALEX2CDE = np.append(data_ALEX2CDE, arrALEX2CDE)
        data_DTGR_TR0 = np.append(data_DTGR_TR0, arrDTGR_TR0)

        BN = BN + numB # increase burst number

        time.sleep(0.05)
        bar()

print('...Burst analysis done!')
print('=====================================================')

# generate dictionary from data arrays
dataOUT_dict = {
    'BN': data_BN,
    'PosT': data_PosT,
    'BIN_T': data_BIN_T,
    'ID': data_ID,
    'IA': data_IA,
    'IA0': data_IA0,
    'TauD': data_TauD,
    'TauA0': data_TauA0,
    'FRET2CDE': data_FRET2CDE,
    'ALEX2CDE': data_ALEX2CDE,
    'DTGR_TR0': data_DTGR_TR0,
}

settings['Algorithm'] = ALGORITHM
settings['Bin_T'] = BIN_T
settings['Threshold'] = THRE_B
settings['Mean_IRF_Donor'] = MEAN_IRF_DONOR
settings['Mean_IRF_Acceptor'] = MEAN_IRF_ACCEPTOR

# save parameters in json file
if not os.path.exists("results"):
    os.makedirs("results")

results_path = os.path.join("results", f"Results_{FOLDER}.json")
settings_path = os.path.join("results", f"Settings_{FOLDER}.json")

with open(results_path, 'w') as f:
    json.dump(dataOUT_dict, f, indent=4, cls=NumpyEncoder)

with open(settings_path, 'w') as f:
    json.dump(settings, f, indent=4)

print('Data saved in result folder!')
