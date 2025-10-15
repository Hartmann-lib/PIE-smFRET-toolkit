# PIE-smFRET-toolkit
Toolkit for single-molecule FRET (smFRET) experiments using Time-Correlated Single Photon Counting (TCSPC) data in PTU format (T3 data, PicoQuant). The scripts support pulsed-interleaved excitation (PIE) of two wavelengths and two detection channels, i.e., donor and acceptor excitations and corresponding detection channels. The measurement folder is supposed to contain multiple PTU files with only few minutes length, preferably 1Min length.

Required Python packages:

- numpy 
- matplotlib
- pandas
- openpyxl
- alive_progress

Usage:

Use the scripts according to there name in the order from S0 to S3. User input like, e.g., channel assignments or filter settings, is required in the parameter section in the beginning of each script.

S0: Defines the microtime windows of donor excitation donor/acceptor emission and acceptor excitation acceptor emission.
In order to load example data for window selection, a folder has to be specified in the parameter section. Furthermore, the number of the donor and acceptor channels have to be set.

S1: Extracts the background counts from the measurement files:
              BD ... average background intensity of the donor channel after donor excitation
              BA ... average background intensity of the acceptor channel after donor excitation
              BA0 ... average background intensity of the acceptor channel after acceptor excitation
Therefore, the measurement folder has to be specified and the settings file derived from script S0 has to be selected. Further settings are:

              BIN_T ... the bin time to calculate the time trace from the photon arrival times
              FRAC_T ... time trace part, which will be shown for background selection
              NUM_F ... number of files to be analyzed
              REG_F ... number of time trace regions per file for background estimation

S2: Analyses the PTU files
