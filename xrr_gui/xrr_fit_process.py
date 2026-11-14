"""The fit runs in its OWN PROCESS so the GUI never shares an interpreter lock with the number crunching (a thread made the window stall for seconds and slowed the fit about 3x).
This module is the child's entry point: it loads the engine, runs analyze_stack, and streams small messages to the GUI through a multiprocessing queue:
    (idx, "log", text)                 engine printout            (idx = the run number: several runs with different seeds can share one queue)
    (idx, "progress", dict)            every 0.25 s: label, stage, step, redchi, elapsed, x / model of the current best, new history entries, new fit starts
    (idx, "done", out) | (idx, "stopped", None) | (idx, "error", traceback)
It imports no GUI code (no tkinter)."""
import sys, io, time, threading, traceback, contextlib
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np


class _MPStop:
    """Looks like a threading.Event to the hook, backed by a multiprocessing.Event the GUI sets."""
    def __init__(self, ev): self.ev = ev
    def is_set(self): return self.ev.is_set()
    def clear(self): pass
    def set(self): self.ev.set()


class _Writer(io.TextIOBase):
    def __init__(self, q, idx): self.q = q; self.idx = idx
    def write(self, s):
        if s:
            self.q.put((self.idx, "log", s))
        return len(s)


def run(payload, ev_q, stop_ev):
    idx = payload.get("idx", 0)
    try:
        import xrr_engine_loader as EL
        hook = EL.FitHook(); hook.stop = _MPStop(stop_ev)
        ns = EL.load_engine(hook, payload.get("notebook"))
        EL.set_options(ns, payload["mode"], payload["fit_range"], payload["seed"], payload["bd"]); ns["_BLANK_NAMES"] = set(payload["blank"])
        ns["CHECK_PROFILE"] = payload.get("check_profile", "calibrated"); ns["SYS_FLOOR"] = float(payload.get("sys_floor", 0.05))      # quality-gate profile (see the notebook notes)
        sent ={"hist": 0, "starts": 0}; done = threading.Event()

        def emit():
            st = dict(hook.state); h = list(hook.history); fs = list(hook.fit_starts); msg = dict(label=hook.label, fit_no=hook.fit_no, hist=h[sent["hist"]:], starts=fs[sent["starts"]:], t=time.time() - hook.t0)
            sent["hist"], sent["starts"] = len(h), len(fs)
            i = st.get("info")
            if st:
                msg.update(stage=st["stage"], step=st["step"], redchi=st["redchi"])
                if i is not None and st.get("u") is not None:
                    try:
                        msg["x"] = np.asarray(i["x"]); msg["model"] = ns["model_counts"](i["full"](st["u"]), i["names"], i["x"], i["lam"], ns["SUBSTRATE_BETA_OVER_DELTA"], i["footprint"])
                    except Exception:
                        pass
            ev_q.put((idx, "progress", msg))

        def sender():
            while not done.is_set():
                time.sleep(0.25)
                if hook.state:
                    emit()
        th = threading.Thread(target=sender, daemon=True); th.start()
        try:
            with contextlib.redirect_stdout(_Writer(ev_q, idx)):
                out = ns["analyze_stack"](payload["scan"], payload["layers"], payload["mats"], payload["sub"], verbose=2, film_id=payload["fid"], compare=payload["compare"])
            done.set(); th.join(1.0); emit(); ev_q.put((idx, "done", out))
        except EL.StopFit:
            done.set(); ev_q.put((idx, "stopped", None))
    except Exception:
        ev_q.put((idx, "error", traceback.format_exc()))
