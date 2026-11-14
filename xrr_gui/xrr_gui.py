"""XRR fitting GUI (prototype 3).  Opens a Rigaku .rasx scan, takes a layer table (add / remove / reorder, substrate flag, blank density = free SLD), fits it with the validated notebook
engine IN SEPARATE PROCESSES (optionally several random seeds in parallel, to show whether the answer depends on the seed), shows the data with the fit progression and the convergence live,
and finishes with the film viewer's XRR panel.  Needs Python 3.9+, numpy, scipy, pandas, matplotlib, Pillow (all installed with matplotlib) and tkinter (included with Python).
usage:  python xrr_gui.py        |   python xrr_gui.py --selftest scan.rasx [--preset NAME] [--mode quick] [--seeds 4]"""
import sys, io, os, queue, time, traceback, argparse, multiprocessing as mp
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import xrr_fit_process
# On Windows every fit process is started by re-importing THIS file as '__mp_main__'. A fit process needs no window, so keep that import light (no Tk, matplotlib, Pillow):
# with several seeds, loading the GUI libraries in every child made the start slow.
IS_CHILD = mp.current_process().name != "MainProcess"
if not IS_CHILD:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from PIL import Image, ImageTk
    import xrr_style, xrr_layers as XL, xrr_engine_loader as EL, xrr_panel, xrr_seeds
    from xrr_layers import Row
    FONT_FAMILY = xrr_style.apply()
    SEED_COLORS = matplotlib.colormaps["tab10"].colors
MAXPAR = max(1, min(8, (os.cpu_count() or 4) - 2))          # seeds run in parallel, leaving two cores for the rest of the machine
COLS = [("role", "", 84), ("name", "name", 90), ("formula", "formula", 80), ("density", "density", 62), ("t_min", "t min", 48), ("t_max", "t max", 48), ("r_lo", "range lo", 58), ("r_hi", "range hi", 58)]


def num(s):
    s = str(s).strip()
    return None if s == "" else float(s)


def stage_key(s):
    return "multi-start" if s.startswith("multi") else ("final polish" if s == "polish" else s)


class App:
    def __init__(self, root):
        self.root = root; root.title("XRR fitting tool"); root.geometry("1500x900"); root.protocol("WM_DELETE_WINDOW", self.close)
        self.ev_q = None; self.ctx = None; self.rs = {}; self.pending = []; self.stop_t = None; self.base = None; self.seed_df = None; self.seed_lines = []
        self.ns = None; self.scan = None; self.scan_path = None; self.out = None; self.panel_fig = None; self._photo = None; self._img = None
        self.rows = [Row(**vars(r)) for r in XL.PRESETS["Fe / Fe2O3 on SrTiO3 (expected)"]]; self.sel = None; self.finished = False; self.running = False
        self.dirty = False; self.last_draw = 0.0; self.sep = []
        self.status = tk.StringVar(value="loading the fitting engine ..."); self._build(); root.update_idletasks()
        self.ns = EL.load_engine(EL.FitHook()); self.status.set("ready: open a .rasx scan"); self._refresh_tree(); self._tick()

    # ------------------------------------------------------------------ layout
    def _build(self):
        main = ttk.Panedwindow(self.root, orient="horizontal"); main.pack(fill="both", expand=True)
        left = ttk.Frame(main, padding=6); right = ttk.Frame(main); main.add(left, weight=0); main.add(right, weight=1)
        f = ttk.LabelFrame(left, text="1. Scan", padding=6); f.pack(fill="x")
        ttk.Button(f, text="Open .rasx ...", command=self.open_scan).grid(row=0, column=0, sticky="w")
        self.scan_lbl = ttk.Label(f, text="no scan loaded", wraplength=420); self.scan_lbl.grid(row=0, column=1, sticky="w", padx=8)
        ttk.Label(f, text="film id").grid(row=1, column=0, sticky="w", pady=(4, 0)); self.fid = tk.StringVar(); ttk.Entry(f, textvariable=self.fid, width=34).grid(row=1, column=1, sticky="w", padx=8, pady=(4, 0))

        f = ttk.LabelFrame(left, text="2. Sample, top (surface) to bottom (substrate)", padding=6); f.pack(fill="x", pady=6)
        self.tree = ttk.Treeview(f, columns=[c[0] for c in COLS], show="headings", height=7, selectmode="browse")
        for k, h, w in COLS:
            self.tree.heading(k, text=h); self.tree.column(k, width=w, anchor="center" if k not in ("name", "formula") else "w")
        self.tree.grid(row=0, column=0, columnspan=6, sticky="ew"); self.tree.bind("<<TreeviewSelect>>", self._on_select)
        self.tree.tag_configure("sub", foreground="#1b4f9c", font=("TkDefaultFont", 9, "bold")); self.tree.tag_configure("ign", foreground="#a0a0a0")
        for i, (t, cmd) in enumerate([("Add layer", self.add_row), ("Duplicate", self.dup_row), ("Remove", self.del_row), ("Up", lambda: self.move(-1)), ("Down", lambda: self.move(1)), ("Set as substrate", self.set_sub)]):
            ttk.Button(f, text=t, command=cmd, width=13 if i == 5 else 9).grid(row=1, column=i, pady=4, padx=1, sticky="ew")
        pf = ttk.Frame(f); pf.grid(row=2, column=0, columnspan=6, sticky="ew")
        self.preset = tk.StringVar(value=list(XL.PRESETS)[0]); ttk.Combobox(pf, textvariable=self.preset, values=list(XL.PRESETS), state="readonly", width=44).pack(side="left")
        ttk.Button(pf, text="Load preset", command=self.load_preset).pack(side="left", padx=4)
        ef = ttk.LabelFrame(f, text="selected row", padding=4); ef.grid(row=3, column=0, columnspan=6, sticky="ew", pady=(6, 0))
        self.v = {k: tk.StringVar() for k in ("name", "formula", "density", "t_min", "t_max", "r_lo", "r_hi")}
        labels = [("name", "name"), ("formula", "formula"), ("density", "density g/cm3 (blank = free)"), ("t_min", "thickness min (A)"), ("t_max", "thickness max (A)"), ("r_lo", "range lo"), ("r_hi", "range hi")]
        for i, (k, t) in enumerate(labels):
            r, c = divmod(i, 2); ttk.Label(ef, text=t).grid(row=r, column=2 * c, sticky="w", padx=2); e = ttk.Entry(ef, textvariable=self.v[k], width=14); e.grid(row=r, column=2 * c + 1, padx=4, pady=1); e.bind("<Return>", lambda _e: self.apply_row())
        ttk.Button(ef, text="Apply to row", command=self.apply_row).grid(row=4, column=0, columnspan=2, pady=4, sticky="w")
        self.hint = ttk.Label(ef, text="", wraplength=420, foreground="#444"); self.hint.grid(row=5, column=0, columnspan=4, sticky="w")
        ttk.Label(f, text="Thickness bounds blank = automatic (from the fringes). Range: density given -> fraction of bulk (default 0.5 - 1.4); density blank -> absolute SLD in 1e-5 A^-2 (default 1 - 10).",
                  wraplength=430, foreground="#555").grid(row=4, column=0, columnspan=6, sticky="w", pady=(4, 0))

        f = ttk.LabelFrame(left, text="3. Fit options", padding=6); f.pack(fill="x")
        self.mode = tk.StringVar(value="quick"); self.lo = tk.StringVar(value="0.2"); self.hi = tk.StringVar(value=""); self.seed = tk.StringVar(value="0"); self.nseeds = tk.StringVar(value="1"); self.bd = tk.StringVar(value="0.05"); self.compare = tk.BooleanVar(value=False)
        ttk.Label(f, text="mode").grid(row=0, column=0, sticky="w"); ttk.Combobox(f, textvariable=self.mode, values=["quick", "standard", "thorough"], state="readonly", width=10).grid(row=0, column=1, padx=4)
        ttk.Label(f, text="2theta from").grid(row=0, column=2, sticky="e"); ttk.Entry(f, textvariable=self.lo, width=6).grid(row=0, column=3)
        ttk.Label(f, text="to (blank = auto)").grid(row=0, column=4, sticky="e"); ttk.Entry(f, textvariable=self.hi, width=6).grid(row=0, column=5)
        ttk.Label(f, text="first seed").grid(row=1, column=0, sticky="w"); ttk.Entry(f, textvariable=self.seed, width=6).grid(row=1, column=1, sticky="w", padx=4)
        ttk.Label(f, text="number of seeds").grid(row=1, column=2, sticky="e"); ttk.Entry(f, textvariable=self.nseeds, width=6).grid(row=1, column=3)
        ttk.Label(f, text="substrate beta/delta").grid(row=1, column=4, sticky="e"); ttk.Entry(f, textvariable=self.bd, width=6).grid(row=1, column=5)
        ttk.Checkbutton(f, text="also refit without each layer (is every layer needed?)", variable=self.compare).grid(row=2, column=0, columnspan=6, sticky="w", pady=(4, 0))
        self.profile = tk.StringVar(value="calibrated"); self.floor = tk.StringVar(value="5")
        ttk.Label(f, text="quality checks").grid(row=3, column=0, sticky="w", pady=(4, 0)); ttk.Combobox(f, textvariable=self.profile, values=["calibrated", "strict"], state="readonly", width=10).grid(row=3, column=1, padx=4, pady=(4, 0))
        ttk.Label(f, text="systematic floor (%)").grid(row=3, column=2, columnspan=2, sticky="e", pady=(4, 0)); ttk.Entry(f, textvariable=self.floor, width=6).grid(row=3, column=4, sticky="w", pady=(4, 0))
        ttk.Label(f, text=f"Several seeds run in parallel (up to {MAXPAR} at once) and show whether the answer depends on the random start of the global search. 'calibrated' judges the fit against counting noise plus the systematic floor; 'strict' uses the original noise-only gates.", wraplength=430, foreground="#555").grid(row=4, column=0, columnspan=6, sticky="w")

        f = ttk.Frame(left, padding=(0, 8)); f.pack(fill="x")
        self.run_btn = ttk.Button(f, text="Run fit", command=self.start_fit); self.run_btn.pack(side="left"); self.stop_btn = ttk.Button(f, text="Stop", command=self.stop_fit, state="disabled"); self.stop_btn.pack(side="left", padx=6)
        ttk.Button(f, text="Save results ...", command=self.save_results).pack(side="left", padx=6)
        ttk.Label(left, textvariable=self.status, wraplength=470, foreground="#1b4f9c").pack(fill="x")

        nb = ttk.Notebook(right); nb.pack(fill="both", expand=True); self.nb = nb
        t1 = ttk.Frame(nb); nb.add(t1, text="Live fit"); t2 = ttk.Frame(nb); nb.add(t2, text="Result panel"); t4 = ttk.Frame(nb); nb.add(t4, text="Seeds"); t3 = ttk.Frame(nb); nb.add(t3, text="Log")
        self.fig = Figure(figsize=(8, 7)); gs = self.fig.add_gridspec(2, 1, height_ratios=[3, 2], hspace=0.34)
        self.ax = self.fig.add_subplot(gs[0]); self.axc = self.fig.add_subplot(gs[1])
        self.ax.set_xlabel(r"$2\theta/\omega$ (deg)"); self.ax.set_ylabel("counts"); self.ax.set_yscale("log"); self.axc.set_xlabel("optimiser step (one point per generation / improvement)"); self.axc.set_ylabel(r"reduced $\chi^2$ of the fit being run"); self.axc.set_yscale("log")
        self.canvas = FigureCanvasTkAgg(self.fig, master=t1); self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.data_out, = self.ax.plot([], [], ".", ms=2, color="0.78", label="data outside the fit range"); self.data_in, = self.ax.plot([], [], ".", ms=2, color="0.35", label="data in the fit range")
        self.model_line, = self.ax.plot([], [], color="crimson", lw=1.4, label="best model so far"); self.ax.legend(loc="lower left", fontsize=8)
        self.conv_lines = {k: self.axc.plot([], [], "o", ms=2.2, color=c, label=k)[0] for k, c in (("global search", "tab:blue"), ("multi-start", "0.55"), ("refine", "darkorange"), ("final polish", "crimson"))}
        self.axc.legend(loc="upper right", fontsize=7, ncol=4)
        self.rbar = ttk.Frame(t2); self.rbar.pack(fill="x")
        ttk.Button(self.rbar, text="Save panel PNG ...", command=lambda: self.save_panel("png")).pack(side="left", padx=4, pady=3); ttk.Button(self.rbar, text="Save panel PDF ...", command=lambda: self.save_panel("pdf")).pack(side="left", pady=3)
        self.rlbl = ttk.Label(t2, text="run a fit to see the result panel", anchor="center"); self.rlbl.pack(fill="both", expand=True); self.rlbl.bind("<Configure>", lambda e: self._show_panel())
        self.seed_txt = tk.Text(t4, height=7, wrap="word", font=("TkDefaultFont", 10)); self.seed_txt.pack(fill="x"); self.seed_txt.insert("end", "Set 'number of seeds' above 1 and run: each seed is an independent fit; this tab then shows the seed-to-seed spread.")
        self.seed_tree = ttk.Treeview(t4, show="headings", height=14); self.seed_tree.pack(fill="both", expand=True, padx=2, pady=4)
        self.log = tk.Text(t3, wrap="none", font=("Consolas", 9)); sb = ttk.Scrollbar(t3, command=self.log.yview); self.log.configure(yscrollcommand=sb.set); sb.pack(side="right", fill="y"); self.log.pack(fill="both", expand=True)

    # ------------------------------------------------------------------ table
    def _refresh_tree(self, select=None):
        self.tree.delete(*self.tree.get_children()); sub = [i for i, r in enumerate(self.rows) if r.substrate]; s = sub[0] if sub else len(self.rows)
        f = lambda x: "" if x is None else f"{x:g}"
        for i, r in enumerate(self.rows):
            tag = "sub" if r.substrate else ("ign" if i > s else "")
            self.tree.insert("", "end", iid=str(i), values=("SUBSTRATE" if r.substrate else ("ignored" if i > s else f"layer {i + 1}"), r.name, r.formula, "free SLD" if r.density is None and not r.substrate else f(r.density), f(r.t_min), f(r.t_max), f(r.r_lo), f(r.r_hi)), tags=(tag,) if tag else ())
        if select is not None and 0 <= select < len(self.rows): self.tree.selection_set(str(select)); self.sel = select
        else: self.sel = None

    def _on_select(self, _e=None):
        s = self.tree.selection()
        if not s: return
        self.sel = int(s[0]); r = self.rows[self.sel]; g = lambda x: "" if x is None else f"{x:g}"
        for k, val in (("name", r.name), ("formula", r.formula), ("density", g(r.density)), ("t_min", g(r.t_min)), ("t_max", g(r.t_max)), ("r_lo", g(r.r_lo)), ("r_hi", g(r.r_hi))): self.v[k].set(val)
        self.hint.config(text=("Blank density: the SLD of this layer is fitted freely between range lo and hi (absolute, 1e-5 A^-2); no formula needed." if r.density is None and not r.substrate else
                               "Density given: it is the bulk reference; the fitted electron density may vary by range lo - hi times bulk." if not r.substrate else "Substrate: needs formula and density; its SLD may vary by +-15 % around bulk in the fit."))

    def apply_row(self):
        if self.sel is None: return
        try:
            r = self.rows[self.sel]; r.name = self.v["name"].get().strip(); r.formula = self.v["formula"].get().strip()
            r.density, r.t_min, r.t_max, r.r_lo, r.r_hi = (num(self.v[k].get()) for k in ("density", "t_min", "t_max", "r_lo", "r_hi"))
        except ValueError:
            messagebox.showerror("Row", "numbers only in the density, thickness and range boxes (or leave blank)"); return
        self._refresh_tree(self.sel); self._on_select()

    def add_row(self):
        i = self.sel if self.sel is not None else next((k for k, r in enumerate(self.rows) if r.substrate), len(self.rows))
        self.rows.insert(i, Row(name=f"layer{len(self.rows)}", formula="", density=None)); self._refresh_tree(i); self._on_select()

    def dup_row(self):
        if self.sel is None: return
        r = Row(**vars(self.rows[self.sel])); r.substrate = False; r.name += "_2"; self.rows.insert(self.sel, r); self._refresh_tree(self.sel); self._on_select()

    def del_row(self):
        if self.sel is None or len(self.rows) <= 1: return
        del self.rows[self.sel]; self._refresh_tree(min(self.sel, len(self.rows) - 1)); self._on_select()

    def move(self, d):
        if self.sel is None: return
        j = self.sel + d
        if 0 <= j < len(self.rows): self.rows[self.sel], self.rows[j] = self.rows[j], self.rows[self.sel]; self._refresh_tree(j); self._on_select()

    def set_sub(self):
        if self.sel is None: return
        for i, r in enumerate(self.rows): r.substrate = (i == self.sel)
        self._refresh_tree(self.sel); self._on_select()

    def load_preset(self):
        self.rows = [Row(**vars(r)) for r in XL.PRESETS[self.preset.get()]]; self._refresh_tree(); self.status.set(f"preset loaded: {XL.describe(self.rows)}")

    # ------------------------------------------------------------------ scan
    def open_scan(self, path=None):
        path = path or filedialog.askopenfilename(title="Rigaku XRR scan", filetypes=[("Rigaku .rasx", "*.rasx"), ("all files", "*.*")])
        if not path: return
        try:
            scan = self.ns["load_rasx"](Path(path), False); self.ns["ensure_reflectivity"](scan)
        except Exception as e:
            messagebox.showerror("Scan", str(e)); return
        self.scan, self.scan_path = scan, Path(path); tt = scan["two_theta"]; notes = self.ns["check_scan"](scan)
        self.scan_lbl.config(text=f"{Path(path).name}\n{len(tt)} points, 2theta {tt.min():.2f}-{tt.max():.2f} deg, step {scan['step']}, {scan['seconds_per_point']:.2f} s/point, max {scan['counts'].max():.3g} counts" + ("".join("\nnote: " + m for m in notes)))
        self.fid.set(self.fid.get() or Path(path).stem.replace(" ", "_")); self._reset_live(); self.status.set("scan loaded; define the sample and press Run fit")

    def _reset_live(self):
        tt, I = self.scan["two_theta"], np.clip(self.scan["counts"], 0.5, None)
        self.data_out.set_data(tt, I); self.data_in.set_data([], []); self.model_line.set_data([], []); self.ax.set_title("")
        for ln in self.conv_lines.values(): ln.set_data([], [])
        for ln in self.seed_lines: ln.remove()
        for s in self.sep: s.remove()
        self.seed_lines = []; self.sep = []; self.rs = {}; self.ax.relim(); self.ax.autoscale_view(); self.axc.set_xlim(0, 10); self.canvas.draw_idle()

    # ------------------------------------------------------------------ fit runs (one process per seed)
    def start_fit(self):
        if self.running: return
        if self.scan is None: messagebox.showinfo("Run", "open a .rasx scan first"); return
        try:
            layers, mats, sub = XL.to_engine(self.rows, self.ns)
            lo, hi = float(self.lo.get()), num(self.hi.get()); seed0 = int(self.seed.get()); bd = float(self.bd.get()); n = int(self.nseeds.get()); floor = float(self.floor.get()) / 100.0
            if not 1 <= n <= 64: raise ValueError("number of seeds must be between 1 and 64")
            if not 0.0 <= floor <= 0.5: raise ValueError("systematic floor must be between 0 and 50 %")
        except Exception as e:
            messagebox.showerror("Cannot run", str(e)); return
        fid = "".join(c if c.isalnum() or c in "._-" else "_" for c in (self.fid.get() or "scan"))
        self.base = dict(scan=self.scan, layers=layers, mats=mats, sub=sub, fid=fid, compare=self.compare.get(), mode=self.mode.get(), fit_range=(lo, hi), bd=bd, check_profile=self.profile.get(), sys_floor=floor,
                         blank=sorted(r.name for r in self.rows if r.density is None and not r.substrate), notebook=self.ns["_NOTEBOOK"])
        self.ns["_BLANK_NAMES"] = set(self.base["blank"]); self.ns["SUBSTRATE_BETA_OVER_DELTA"] = bd           # the panel is drawn in this process
        self._reset_live(); self.out = None; self.seed_df = None; self.finished = False; self.running = True; self.stopping = False; self.stop_t = None
        self.run_btn.config(state="disabled"); self.stop_btn.config(state="normal"); self.log.delete("1.0", "end"); self.t_start = time.time()
        self.rs = {i: dict(seed=seed0 + i, hist=[], starts=[], live=None, x=None, model=None, main_chi=np.inf, state="queued", proc=None, stop_ev=None, out=None, err=None) for i in range(n)}
        self.ctx = mp.get_context("spawn"); self.ev_q = self.ctx.Queue(); self.pending = list(range(n)); self._launch()
        self.status.set(f"starting {n} fit process(es): {XL.describe(self.rows)} ({self.mode.get()} mode) ...")

    def _launch(self):
        if self.pending and sum(1 for r in self.rs.values() if r["state"] == "running") < MAXPAR:          # one process per call (_tick calls this every 150 ms): starting them all at once froze the window for seconds
            i = self.pending.pop(0); r = self.rs[i]; r["stop_ev"] = self.ctx.Event(); r["state"] = "running"
            r["proc"] = self.ctx.Process(target=xrr_fit_process.run, args=(dict(self.base, seed=r["seed"], idx=i), self.ev_q, r["stop_ev"]), daemon=True); r["proc"].start()

    def stop_fit(self):
        if self.running:
            self.stopping = True; self.pending = []; self.stop_t = time.time(); self.status.set("stopping ...")
            for r in self.rs.values():
                if r["state"] == "running": r["stop_ev"].set()

    def close(self):
        for r in self.rs.values():
            try:
                if r["proc"] is not None and r["proc"].is_alive(): r["proc"].terminate()
            except Exception: pass
        self.root.destroy()

    def _on_progress(self, idx, m):
        r = self.rs[idx]; r["hist"] += list(m["hist"]); r["starts"] += list(m["starts"]); r["live"] = m
        if "x" in m and "model" in m: r["x"], r["model"] = m["x"], m["model"]
        if m.get("label") == "main fit" and "redchi" in m: r["main_chi"] = min(r["main_chi"], m["redchi"])
        self.dirty = True

    def _draw_live(self):
        runs = list(self.rs.values()); n = len(runs); shown = [r for r in runs if r["x"] is not None]
        if not shown: return
        best = min(shown, key=lambda r: r["main_chi"]); x = best["x"]; m = best["live"]
        self.model_line.set_data(x, np.clip(best["model"], 0.5, None)); tt, I = self.scan["two_theta"], np.clip(self.scan["counts"], 0.5, None); inn = (tt >= x.min() - 1e-9) & (tt <= x.max() + 1e-9)
        self.data_in.set_data(tt[inn], I[inn]); self.data_out.set_data(tt[~inn], I[~inn])
        ndone = sum(r["state"] in ("done", "stopped", "error") for r in runs); nrun = sum(r["state"] == "running" for r in runs)
        if n == 1:
            self.ax.set_title(f"fit {m['fit_no']}: {m['label']}   (reduced $\\chi^2$ {m.get('redchi', float('nan')):.4g}, {m['t']:.0f} s)", fontsize=9)
            h = best["hist"]
            if h:
                ys = np.array([a[3] for a in h]); st = np.array([stage_key(a[1]) for a in h]); xs = np.arange(len(h))
                for k, ln in self.conv_lines.items(): ln.set_data(xs[st == k], ys[st == k])
                self.axc.set_xlim(-0.02 * len(h), len(h) * 1.03 + 1); self.axc.set_ylim(max(ys.min() * 0.7, 1e-3), ys.max() * 1.6)
                for s in self.sep: s.remove()
                self.sep = []
                for idx, lab in best["starts"]:
                    self.sep.append(self.axc.axvline(idx, color="0.8", lw=0.7, ls="--", zorder=0))
                    self.sep.append(self.axc.text(idx, 1.0, " " + lab.split(":")[0].replace("systematic ", "").replace("stack without a layer", "w/o layer")[:16], transform=self.axc.get_xaxis_transform(), fontsize=6, color="0.4", rotation=90, va="top", ha="left"))
            self.status.set(f"fit {m['fit_no']} ({m['label']}), {m.get('stage', '')}: reduced chi2 {m.get('redchi', float('nan')):.4g}; {m['t']:.0f} s")
        else:
            bi = [i for i, r in self.rs.items() if r is best][0]
            self.ax.set_title(f"best of {n} seeds so far: seed {best['seed']} (main-fit reduced $\\chi^2$ {best['main_chi']:.4g}); {ndone} finished, {nrun} running", fontsize=9)
            while len(self.seed_lines) < n: self.seed_lines.append(self.axc.plot([], [], lw=0.9)[0])
            ymin, ymax, xmax = np.inf, 0.0, 1
            for i, r in self.rs.items():
                h = r["hist"]; cut = r["starts"][1][0] if len(r["starts"]) > 1 else len(h)                   # the main fit only (later fits are refits)
                h = h[:cut]
                if not h: continue
                ys = np.array([a[3] for a in h]); ln = self.seed_lines[i]; ln.set_data(np.arange(len(h)), ys); ln.set_color(SEED_COLORS[i % 10]); ln.set_label(f"seed {r['seed']}")
                ymin, ymax, xmax = min(ymin, ys.min()), max(ymax, ys.max()), max(xmax, len(h))
            if ymax > 0: self.axc.set_xlim(0, xmax * 1.03 + 1); self.axc.set_ylim(max(ymin * 0.7, 1e-3), ymax * 1.6)
            if n <= 10: self.axc.legend(handles=self.seed_lines[:n], loc="upper right", fontsize=7, ncol=min(n, 5))
            self.status.set(f"seeds: {ndone}/{n} finished, {nrun} running; best main-fit reduced chi2 {best['main_chi']:.4g} (seed {best['seed']}); {time.time() - self.t_start:.0f} s")
        self.ax.relim(); self.ax.autoscale_view(); self.canvas.draw_idle()

    def _tick(self):
        if self.ev_q is not None:
            try:
                for _ in range(600):
                    idx, kind, data = self.ev_q.get_nowait()
                    if kind == "log":
                        if len(self.rs) > 1 and data.strip(): data = f"[seed {self.rs[idx]['seed']}] " + data
                        self.log.insert("end", data); self.log.see("end")
                    elif kind == "progress": self._on_progress(idx, data)
                    else: self._run_end(idx, kind, data)
            except queue.Empty: pass
            except (EOFError, OSError): pass
        if self.running and self.pending and not self.stopping: self._launch()
        if self.running and self.dirty and time.time() - self.last_draw > 0.5:
            self._draw_live(); self.dirty = False; self.last_draw = time.time()
        if self.running and self.stop_t is not None and time.time() - self.stop_t > 5:                    # cooperative stop did not answer: end the processes
            for i, r in self.rs.items():
                if r["state"] == "running" and r["proc"] is not None and r["proc"].is_alive(): r["proc"].terminate(); self._run_end(i, "stopped", None)
        if self.running:
            for i, r in self.rs.items():                                                                 # a process that died without a message
                if r["state"] == "running" and r["proc"] is not None and not r["proc"].is_alive():
                    r.setdefault("_dead", time.time())
                    if time.time() - r["_dead"] > 2.0: self._run_end(i, "error", "the fit process ended unexpectedly")
        self.root.after(150, self._tick)

    def _run_end(self, idx, kind, data):
        r = self.rs[idx]
        if r["state"] != "running": return
        r["state"] = {"done": "done", "stopped": "stopped"}.get(kind, "error")
        if kind == "done": r["out"] = data
        elif kind == "error": r["err"] = str(data); self.log.insert("end", f"[seed {r['seed']}] ERROR\n{data}\n")
        self._launch()
        if not self.pending and all(x["state"] != "running" for x in self.rs.values()): self._finalize()

    def _finalize(self):
        self.running = False; self.run_btn.config(state="normal"); self.stop_btn.config(state="disabled"); self.finished = True
        if self.dirty: self._draw_live()
        outs = {r["seed"]: r["out"] for r in self.rs.values() if r["out"] is not None}; n = len(self.rs)
        if not outs:
            errs = [r["err"] for r in self.rs.values() if r["err"]]
            if errs: messagebox.showerror("Fit failed", errs[0].strip().split("\n")[-1]); self.status.set("fit failed (see the Log tab)")
            else: self.status.set("stopped by the user")
            return
        extra = []
        if len(outs) > 1:
            self.seed_df = xrr_seeds.seed_table(outs, xrr_panel.fix_checks); stable, lines, _d = xrr_seeds.stability(self.seed_df)
            best_seed = xrr_seeds.solution_groups(self.seed_df)[0]["seeds"][0]; self.out = outs[best_seed]; self._show_seeds(lines, best_seed)
            extra = [f"lowest-chi-square solution shown (seed {best_seed})"] + lines
            summary = f"{len(outs)} seeds done in {time.time() - self.t_start:.0f} s: " + ("stable across seeds" if stable else "NOT stable across seeds (see the Seeds tab)")
        else:
            self.out = next(iter(outs.values())); summary = f"done in {time.time() - self.t_start:.0f} s: {xrr_panel.fix_checks(self.out['result'])[1]}"
        self.extra_lines = extra; self._show_final_live(); self._show_result(); self.status.set(summary)

    def _show_seeds(self, lines, best_seed):
        df = self.seed_df; self.seed_txt.delete("1.0", "end"); self.seed_txt.insert("end", "\n".join(lines) + f"\n\nThe result panel shows seed {best_seed} (lowest reduced chi-square; seeds that ended in the same solution give the same fit).")
        cols = [c for c in df.columns if not c.endswith(("_SLD_1e-5", "_density_frac_of_bulk"))]; self.seed_tree.config(columns=cols)
        for c in cols: self.seed_tree.heading(c, text=c.replace("_thickness_A", " (A)").replace("_roughness_A", " rough (A)").replace("reduced_chi2", "chi2").replace("total_A", "total (A)")); self.seed_tree.column(c, width=130 if c == "verdict" else 82, anchor="center")
        self.seed_tree.delete(*self.seed_tree.get_children())
        for _, row in df.iterrows(): self.seed_tree.insert("", "end", values=[(f"{row[c]:.4g}" if isinstance(row[c], float) else row[c]) for c in cols])

    def _show_final_live(self):
        """After the run the live plot shows the MAIN fit (the one in the result panel), not whichever refit ran last."""
        m = self.out["main"]; x = np.asarray(m["x"]); mod = self.ns["model_counts"](m["v"], m["names"], x, self.scan["wavelength"], self.ns["SUBSTRATE_BETA_OVER_DELTA"], m["footprint"])
        tt, I = self.scan["two_theta"], np.clip(self.scan["counts"], 0.5, None); inn = (tt >= x.min() - 1e-9) & (tt <= x.max() + 1e-9)
        self.model_line.set_data(x, np.clip(mod, 0.5, None)); self.data_in.set_data(tt[inn], I[inn]); self.data_out.set_data(tt[~inn], I[~inn])
        self.ax.set_title(f"final: main fit, reduced $\\chi^2$ {m['redchi']:.4g} (same curve as the result panel)" + ("; the refits that follow estimate the systematic error" if len(self.rs) == 1 else ""), fontsize=9); self.canvas.draw_idle()

    # ------------------------------------------------------------------ result
    def _show_result(self):
        fig = Figure(figsize=(13.5, 9.5)); xrr_panel.draw_stack_panel(fig, self.out, self.scan, self.ns, title_prefix=self.fid.get() or None, extra_lines=getattr(self, "extra_lines", ())); self.panel_fig = fig
        buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=110); buf.seek(0); self._img = Image.open(buf).convert("RGB"); self._show_panel(); self.nb.select(1 if self.seed_df is None else 2)

    def _show_panel(self):
        if self._img is None: return
        w, h = max(self.rlbl.winfo_width(), 200), max(self.rlbl.winfo_height(), 200); s = min(w / self._img.width, h / self._img.height)
        self._photo = ImageTk.PhotoImage(self._img.resize((max(int(self._img.width * s), 1), max(int(self._img.height * s), 1)), Image.LANCZOS)); self.rlbl.config(image=self._photo, text="")

    def save_panel(self, kind):
        if self.panel_fig is None: messagebox.showinfo("Save", "no result yet"); return
        p = filedialog.asksaveasfilename(defaultextension="." + kind, initialfile=f"{self.fid.get() or 'xrr'}_panel.{kind}")
        if p: self.panel_fig.savefig(p, dpi=200 if kind == "png" else None); self.status.set(f"saved {p}")

    def save_results(self, folder=None):
        if self.out is None: messagebox.showinfo("Save", "no result yet"); return
        folder = folder or filedialog.askdirectory(title="Folder for the result files")
        if not folder: return
        folder = Path(folder); names = self.ns["export_fit_details"](self.out, folder); base = self.out["result"]["base_name"]; self.ns["parameter_table"](self.out).to_csv(folder / f"{base}_parameters.csv", index=False)
        self.panel_fig.savefig(folder / f"{base}_panel.png", dpi=200); self.panel_fig.savefig(folder / f"{base}_panel.pdf"); k = len(names) + 3
        if self.seed_df is not None:
            self.seed_df.to_csv(folder / f"{base}_seeds.csv", index=False); (folder / f"{base}_seeds_summary.txt").write_text("\n".join(xrr_seeds.stability(self.seed_df)[1]), encoding="utf-8"); k += 2
        self.status.set(f"saved {k} files in {folder}")


def selftest(a):
    root = tk.Tk(); app = App(root); t0 = time.time()
    app.open_scan(a.selftest); app.preset.set(a.preset); app.load_preset(); app.mode.set(a.mode); app.compare.set(a.compare); app.nseeds.set(str(a.seeds)); app.profile.set(a.profile); app.floor.set(str(a.floor))
    if a.hi: app.hi.set(a.hi)
    app.start_fit(); shots = 0; gaps = []; bigs = []; last = time.time()
    while not app.finished and time.time() - t0 < a.timeout:
        root.update(); time.sleep(0.01); now = time.time(); gaps.append(now - last); last = now
        if gaps[-1] > 0.4: bigs.append((round(now - t0), round(gaps[-1], 2)))
        if a.shot and app.running and shots < 3 and time.time() - t0 > 8 + 14 * shots: _grab(root, f"{a.shot}_live_{shots}.png"); shots += 1
    root.update()
    print("event-loop gaps above 0.4 s (t, gap):", bigs[:20], "| count", len(bigs), "| 99th percentile %.3f s" % np.percentile(gaps, 99))
    print("finished:", app.finished, "| wall time", round(time.time() - t0), "s | status:", app.status.get())
    if app.seed_df is not None:
        import pandas as pd
        pd.set_option("display.width", 250, "display.max_columns", 30); print(app.seed_df.round(3).to_string(index=False)); print("\n".join(xrr_seeds.stability(app.seed_df)[1]))
    if app.out is not None:
        r = app.out["result"]; print("shown fit: total thickness", round(r["thickness_A"], 2), "| reduced chi2", round(r["reduced_chi2"], 1), "|", "; ".join(f"{L['name']} {L['thickness_A']:.1f} A {L['density_frac_of_bulk']*100:.0f}%" for L in reversed(r["layers"])))
        outdir = Path(a.outdir); outdir.mkdir(exist_ok=True, parents=True); app.save_results(outdir); app.fig.savefig(outdir / "live_fit_tab.png", dpi=110)
        if a.shot: root.update(); _grab(root, a.shot + "_result_tab.png"); app.nb.select(0); root.update(); _grab(root, a.shot + "_final_live_tab.png"); app.nb.select(2); root.update(); _grab(root, a.shot + "_seeds_tab.png")
        print("saved to", outdir)
    app.close()


def _grab(root, path):
    try:
        from PIL import ImageGrab
        root.update(); x, y, w, h = root.winfo_rootx(), root.winfo_rooty(), root.winfo_width(), root.winfo_height(); ImageGrab.grab(bbox=(x, y, x + w, y + h)).save(path); print("screenshot", path)
    except Exception as e:
        print("screenshot failed:", e)


if __name__ == "__main__":
    mp.freeze_support()
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest"); ap.add_argument("--preset", default=list(XL.PRESETS)[0]); ap.add_argument("--mode", default="quick"); ap.add_argument("--compare", action="store_true"); ap.add_argument("--hi", default="")
    ap.add_argument("--seeds", type=int, default=1); ap.add_argument("--profile", default="calibrated"); ap.add_argument("--floor", type=float, default=5.0); ap.add_argument("--timeout", type=float, default=1800); ap.add_argument("--outdir", default="selftest_out"); ap.add_argument("--shot", default="")
    a = ap.parse_args()
    if a.selftest: selftest(a)
    else:
        root = tk.Tk(); App(root); root.mainloop()
