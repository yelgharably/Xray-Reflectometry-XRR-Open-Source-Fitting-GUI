# XRR fitting tool

Fits X-ray reflectivity (XRR) scans from Rigaku `.rasx` files with the Parratt formalism, for a single film or a **stack of layers of different materials** on a substrate. It returns every parameter
with uncertainties and a list of quality checks that say whether the result can be trusted.

Two ways to use it, both running the same fitting engine:

| | |
|---|---|
| `XRR_Fitting_Tool_v2.ipynb` | notebook: edit a settings cell (files, substrate, layer table), run, results go to `xrr_results/` |
| `xrr_gui/` | desktop GUI: open a `.rasx`, build the layer table with buttons, watch the fit and its convergence live, get the result panel; optional multi-seed robustness check |

## Quick start
```
pip install -r requirements.txt
python xrr_gui/xrr_gui.py          # the GUI
jupyter notebook XRR_Fitting_Tool_v2.ipynb
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
Prototype. The notebook's engine was validated on simulated data (single film and stack self-tests, run them once on your machine) and on measured films; the GUI reproduces the notebook's results exactly for the same settings.
Bundled font: Latin Modern Roman (GUST Font License, see `xrr_gui/fonts/`).
