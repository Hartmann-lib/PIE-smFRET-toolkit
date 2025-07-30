# PIE-smFRET-toolkit
Toolkit for single-molecule FRET (smFRET) experiments using Time-Correlated Single Photon Counting (TCSPC) data in PTU format (T3 data, PicoQuant). The scripts support pulsed-interleaved excitation (PIE) of two wavelengths and two detection channels, i.e., donor and acceptor excitations and corresponding detection channels.

import os
import tkinter as tk
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from tkinter import filedialog, Tk, messagebox
import json
from scripts.Scatter2Density import Scatter2Density
import pandas as pd
