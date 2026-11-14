"""Loads the validated fitting engine straight from the XRR fitting notebook (all code cells before '## 7. Run on your files', executed unchanged) and wraps the scipy optimisers and
`fit_stack` so the GUI can watch the fit (which fit, convergence, current best model) and stop it.  The fitting maths is NOT modified: the wrappers only observe.

Fits that one analysis runs (in this order): the MAIN fit of the whole stack; optionally the stack WITHOUT each layer in turn; then the systematic-error REFITS of the main fit (a different fit
range, footprint removed).  Each full fit has the stages global differential evolution (one point per generation), multi-start least squares, final polish; a refit only has one polish."""
import inspect, json, os, tempfile, threading, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
NOTEBOOKS = [HERE / "XRR_Fitting_Tool_v2.ipynb", HERE.parent / "XRR_Fitting_Tool_v2.ipynb", HERE.parent / "XRR_Fitting_Tool.ipynb", HERE / "XRR_Fitting_Tool.ipynb"]
POPMAX = {"quick": (10, 400), "standard": (20, 1500), "thorough": (35, 3000)}


class StopFit(Exception):
    """Raised inside the optimiser when the user presses Stop."""


class FitHook:
    """State of the running fit. Written by the engine's thread, read by whoever forwards it to the GUI (xrr_fit_process sends it across a process boundary)."""

    def __init__(self):
        self.stop = threading.Event()        # anything with .is_set() works (the worker process passes a multiprocessing.Event wrapper)
        self.reset()

    def reset(self):
        try:
            self.stop.clear()
        except Exception:
            pass
        self.history = []            # (fit_no, stage, step, reduced chi2): one entry per DE generation / thinned least-squares improvement
        self.fit_starts = []         # (index into history, label) at the start of every fit
        self.state = {}              # latest: fit_no, stage, step, redchi, u, info
        self.fit_no = 0
        self.label = ""
        self.main_n = None           # number of layers of the main fit
        self.stage = ""
        self._ls_calls = 0
        self._had_de = False
        self._last = 0.0
        self._best = np.inf
        self.t0 = time.time()

    def _push(self, stage, step, redchi, u, info, force=False):
        self.stage = stage
        now = time.time()
        if force or now - self._last > 0.08:
            self.history.append((self.fit_no, stage, step, redchi))
            self._last = now
        self.state = dict(fit_no=self.fit_no, stage=stage, step=step, redchi=redchi, u=None if u is None else np.array(u, float), info=info, t=now - self.t0, label=self.label)


def _info(func):
    """names / x / y / lam / footprint / full() of the fit currently running, read from the optimiser's closure (read-only)."""
    try:
        res = inspect.getclosurevars(func).nonlocals.get("resid", func)
        nl = inspect.getclosurevars(res).nonlocals
        return dict(full=nl["full"], names=nl["names"], x=nl["x"], y=nl["y"], lam=nl["lam"], footprint=nl["footprint"])
    except Exception:
        return None


def load_engine(hook: FitHook, notebook: Path = None):
    path = Path(notebook) if notebook else next((p for p in NOTEBOOKS if p.exists()), None)
    if path is None:
        raise FileNotFoundError("XRR_Fitting_Tool_v2.ipynb (or XRR_Fitting_Tool.ipynb) not found next to the GUI folder")
    nb = json.load(open(path, encoding="utf-8")); cells = nb["cells"]
    stop = next(i for i, c in enumerate(cells) if c["cell_type"] == "markdown" and "".join(c["source"]).lstrip().startswith("## 7. Run"))
    ns = {"__name__": "xrr_engine"}
    old = os.getcwd(); tmp = tempfile.mkdtemp(prefix="xrr_engine_")
    os.chdir(tmp)
    try:
        for c in cells[:stop]:
            if c["cell_type"] == "code":
                exec("".join(c["source"]), ns)
    finally:
        os.chdir(old)
    ns["FILES"] = []; ns["LAYERS"] = None; ns["PRESET"] = None
    orig_de, orig_ls, orig_fs = ns["differential_evolution"], ns["least_squares"], ns["fit_stack"]

    def fit_stack(tt, I, lam, layers, mats, sub, lo, hi, t_bounds, scale_bounds, bg_max, footprint=True, start=None, seed=0):
        hook.fit_no += 1; hook._ls_calls = 0; hook._had_de = False; hook._best = np.inf
        if hook.main_n is None:
            hook.main_n = len(layers); hook.label = "main fit"
        elif start is None:
            hook.label = f"stack without a layer ({len(layers)} of {hook.main_n} layers)"
        else:
            hook.label = f"systematic refit: {lo:.2f}-{hi:.2f} deg" + ("" if footprint else ", footprint removed")
        hook.fit_starts.append((len(hook.history), hook.label))
        return orig_fs(tt, I, lam, layers, mats, sub, lo, hi, t_bounds, scale_bounds, bg_max, footprint=footprint, start=start, seed=seed)

    def de(func, bounds, **kw):
        info = _info(func); dof = max(len(info["x"]) - len(bounds), 1) if info else 1
        hook._had_de = True; gen = {"n": 0}

        def cb(intermediate_result):
            if hook.stop.is_set():
                raise StopFit()
            gen["n"] += 1
            hook._push("global search", gen["n"], float(intermediate_result.fun) / dof, intermediate_result.x, info)
            return False
        kw["callback"] = cb
        return orig_de(func, bounds, **kw)

    def ls(fun, x0, **kw):
        info = _info(fun); hook._ls_calls += 1; calls = {"n": 0}
        multistart = kw.get("max_nfev") is not None
        stage = "multi-start" if multistart else ("polish" if not hook._had_de else ("refine" if hook._ls_calls == 1 else "final polish"))
        dof = max(len(info["x"]) - len(x0), 1) if info else 1

        def f(u):
            if hook.stop.is_set():
                raise StopFit()
            r = fun(u); calls["n"] += 1; c = float(np.dot(r, r)) / dof
            improved = c < hook._best
            if improved:
                hook._best = c
            if improved or (not multistart and calls["n"] % 5 == 0):
                hook._push(stage, calls["n"], hook._best, u if improved else None, info)
            return r
        return orig_ls(f, x0, **kw)

    ns["differential_evolution"], ns["least_squares"], ns["fit_stack"] = de, ls, fit_stack
    ns["_NOTEBOOK"] = str(path)
    return ns


def set_options(ns, mode, fit_range=(0.2, None), seed=0, beta_over_delta=0.05, thickness_bounds=None):
    """The notebook computes POP / MAXIT once at load; set them (and the other module-level options the engine reads) before each run."""
    ns["MODE"] = mode; ns["POP"], ns["MAXIT"] = POPMAX[mode]
    ns["FIT_RANGE_2THETA"] = fit_range; ns["SEED"] = int(seed); ns["SUBSTRATE_BETA_OVER_DELTA"] = float(beta_over_delta); ns["THICKNESS_BOUNDS"] = thickness_bounds
    ns["VERBOSITY"] = 2
