# XRR fitting GUI (prototype 3)

Open a Rigaku `.rasx` XRR scan, describe the sample as a layer table, fit it with the same engine as `XRR_Fitting_Tool.ipynb`, and watch the fit live. The result tab is the FINCH film viewer's XRR panel.

**Run:** `python xrr_gui.py` (Python 3.9+; numpy, scipy, pandas, matplotlib, Pillow; tkinter ships with Python, on some Linux installs `sudo apt install python3-tk`).
The folder next to this one must contain `XRR_Fitting_Tool_v2.ipynb` (or `XRR_Fitting_Tool.ipynb`): the GUI loads the fitting engine from the notebook's code cells unchanged.

**Layer table** (top = surface, bottom = substrate, as in the real sample):
- `Add layer` inserts above the selected row; `Duplicate`, `Remove`, `Up`, `Down` reorder; `Set as substrate` marks the substrate (rows below it are ignored). Edit a row in the form and press `Apply to row`.
- Density given: bulk density in g/cm3; the fitted electron density may vary by the range factors relative to bulk (default 0.5-1.4).
- **Density blank = free SLD**: no formula needed; the layer's SLD is fitted freely between range lo and hi, absolute, in 1e-5 A^-2 (default 1-10; SrTiO3 4.0, Fe 6.2, Fe2O3 4.2).
- Thickness min/max blank = automatic (from the fringes). The substrate needs formula and density (its SLD may vary by +-15 %).
- Presets are examples; a one-row table is the single-film case (the notebook's automatic 'film + upper layer' model choice is not in the GUI).

**Fit:** modes quick (about 1 min), standard, thorough. Every fit runs in its own process, so the window stays responsive. The live tab shows the data with the best model so far (data outside the fit range in light grey) and the reduced chi-square against the optimiser step; the title says which fit is being run (the main fit, the stack without a layer, or a systematic refit), and after the run it switches to the main fit, the same curve as the result panel. `Stop` aborts the fit.
`Also refit without each layer` adds the evidence for whether every layer is needed (slower). `Save results...` writes the notebook's files (details JSON, parameters, quality, report) plus the panel as PNG and PDF (and the seed table when several seeds were run).

**Several seeds (robustness check):** set `number of seeds` above 1. Each seed is an independent fit (different random start of the global search); up to 8 run in parallel. The Seeds tab then lists every seed (chi-square, total thickness, each layer), how many DISTINCT solutions were found, and a stability verdict with fixed limits: total thickness spread < 5 %, every layer thickness spread < 15 %, chi-square within a factor 1.5. The result panel shows the lowest-chi-square solution with these lines added. A model that is not stable across seeds is under-determined by the scan: report the spread, not one fit. Use at least 6 seeds: with only 2 or 3 the seeds can agree by chance (for the Fe/STO scan seeds 0 and 1 agree while other seeds find a second minimum), so a 'stable' verdict from a handful of seeds is weak evidence.

**Quality checks:** the checks that test the parameters (stability across refits and seeds, parameters on a bound, every layer needed, fringe agreement, dynamic range) are always strict. The two fit-quality checks have two profiles (`quality checks` in the fit options): `calibrated` (default) judges the misfit against counting noise PLUS a relative systematic floor (default 5 %, `chi2_sys` below 3 passes) and counts residual structure only when the misfit also exceeds that floor; `strict` is the original noise-only gate (reduced chi-square below 15, lag-1 autocorrelation below 0.6), which fails almost every high-count scan. The raw reduced chi-square is always shown. See the notes at the end of the notebook for the evidence behind the calibration. A failed reliability or fit-quality check is a *concern* and is counted in the verdict (`REVIEW: n concern(s)`); the density-plausibility check is only a *note* (the apparent density depends on the assumed formula and on porosity) and is listed but not counted (`PASS: no concerns, 1 note(s)`).

**Self-test without clicking:** `python xrr_gui.py --selftest scan.rasx --preset "Fe / FeTiO3 interlayer on SrTiO3" --mode quick --outdir out`.

Fonts: Latin Modern Roman (GUST Font License, in `fonts/`) so the panel looks the same on every machine.
Limits: the numbers are only as good as the model; check the quality verdict and failed checks in the panel before quoting a result.
