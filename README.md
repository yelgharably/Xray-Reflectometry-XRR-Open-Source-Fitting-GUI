# XRR fitting tool

Fits X-ray reflectivity (XRR) scans from Rigaku `.rasx` files with the Parratt formalism, for a single film or a **stack of layers of different materials** on a substrate. It returns every parameter
with uncertainties and a list of quality checks that say whether the result can be trusted.

Two ways to use it, both running the same fitting engine:

| | |
|---|---|
| `XRR_Fitting_Tool_v2.ipynb` | notebook: edit a settings cell (files, substrate, layer table), run, results go to `xrr_results/` |
| `xrr_gui/` | desktop GUI: open a `.rasx`, build the layer table with buttons, watch the fit and its convergence live, get the result panel; optional multi-seed robustness check |

## Setup
You need **Python 3.9 or newer** (developed and tested on Python 3.13). Check what you have:
```
python --version
```
If it is older than 3.9, install a current Python from https://www.python.org/downloads/ (on Windows tick "Add python.exe to PATH" in the installer).

Create an isolated environment and install the dependencies (numpy, scipy, pandas, matplotlib, Pillow and Jupyter):

Windows (PowerShell):
```
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```
macOS / Linux:
```
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```
Without a virtual environment, the same packages in one line: `pip install numpy scipy pandas matplotlib Pillow jupyter`.

The GUI also needs tkinter, which ships with Python. Check it with `python -c "import tkinter"`; if that fails on Linux, install it with `sudo apt install python3-tk`.
Tested with numpy 2.4, scipy 1.15, pandas 2.2 and matplotlib 3.10 (the code needs scipy 1.11 or newer).

## Quick start
```
python xrr_gui/xrr_gui.py                     # the GUI
jupyter notebook XRR_Fitting_Tool_v2.ipynb    # or the notebook
```
The GUI loads the fitting engine from the notebook's code cells, so keep `XRR_Fitting_Tool_v2.ipynb` in the folder above `xrr_gui/`. Details of the GUI are in [xrr_gui/README.md](xrr_gui/README.md).

## What the fit does
Global differential evolution, then multi-start and trust-region least squares, on the Poisson-weighted residual of log10(counts). Free parameters per layer: thickness, SLD, absorption, roughness of its top interface;
global: scale, background, 2theta offset, substrate SLD, footprint angle, angular resolution. Uncertainties come from the covariance scaled by the reduced chi-square; a systematic error comes from refits under different
choices (fit range, footprint). Quality checks include parameters on bounds, plausible densities, stability across refits, reduced chi-square and residual structure.

## Read this before quoting a number
* XRR measures electron density, so a layer's "density" is only as meaningful as its assumed formula; two different oxides can be indistinguishable.
* Check the verdict and the failed checks. A model that does not fit (large reduced chi-square, structured residuals, parameters on bounds) is missing something; do not quote its layer parameters.
* Use the **multi-seed option** (GUI): if different random starts give different layer structures or thicknesses, the scan does not determine the model. Report the spread, not a single fit. Use at least 6 seeds; two or three can agree by chance.
* Test scans that were not specular reflectivity are refused rather than fitted.

## Status
Prototype. The notebook's engine was validated (single film and stack self-tests, run them once on your machine) on measured films and compared to the Rigaku SmartLab Studio software; the GUI reproduces the notebook's results exactly for the same settings.
Bundled font: Latin Modern Roman (GUST Font License, see `xrr_gui/fonts/`).

<img width="1492" height="923" alt="image" src="https://github.com/user-attachments/assets/2656b581-0b8e-4846-b951-462845e51f05" />

<img width="927" height="652" alt="image" src="https://github.com/user-attachments/assets/c97d26cc-64ee-4e49-940f-a7eb819e185d" />
