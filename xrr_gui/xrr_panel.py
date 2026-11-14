"""The film viewer's XRR panel for a layer-stack fit, drawn from an in-memory engine result (same layout and drawing code as finch_db.plots._plot_xrr_stack_right /
scripts/plot_xrr_stack_viewer_style.py): data + fit, I q^4, normalised residuals on the left; model equations, fitted parameters (red = at a bound), derived quantities and quality on the right,
and an SLD depth profile inset.  For 1-3 layers the right column is the viewer's (four parameter rows per layer); for 4 or more layers it uses one table row per layer, the summary lines are wrapped,
and the figure becomes taller if the content needs the room.  Call xrr_style.apply() first so the fonts match the viewer.  All LaTeX fragments are raw strings (a bare backslash-t would be a tab)."""
import re
import textwrap
import numpy as np
from scipy.special import erf

EQ = [
    r"$q=\frac{4\pi}{\lambda}\sin\theta,\quad \theta=\frac{2\theta_{\mathrm{meas}}-2\theta_0}{2}$",
    r"$k_j=\sqrt{(q/2)^2-4\pi(\rho_j-\rho_0)},\quad \rho_j=(\mathrm{SLD}_j-i\beta_j)\,10^{-5}\,\mathrm{\AA}^{-2}$",
    r"$r_{j,j+1}=\frac{k_j-k_{j+1}}{k_j+k_{j+1}}\;e^{-2k_jk_{j+1}\sigma_j^{2}}$",
    r"$X_j=\frac{r_{j,j+1}+X_{j+1}\,e^{2ik_{j+1}d_{j+1}}}{1+r_{j,j+1}X_{j+1}\,e^{2ik_{j+1}d_{j+1}}},\quad X_{N}=0,\quad R=|X_0|^{2}$",
    r"$I(2\theta)=\left[\,10^{\log I_0}\,\mathrm{erf}\!\left(\frac{2\theta}{\theta_{fp}}\right)R(q)\,\right]\ast G(\sigma_{res})+I_{bg}$",
]


def fix_checks(R):
    """Post-process the engine's check list for display. (1) The density-plausibility check is advisory (a NOTE, not a concern); results written before that split carry it as 'NOT OK', so convert it.
    (2) That check reports nan for free-density layers (blank density): drop those entries. The verdict text is rebuilt from the statuses: concerns ('NOT OK') are counted, notes are not.
    Returns (checks, verdict_text, note or None)."""
    out, skipped = [], []
    for c in R.get("check_details", []):
        c = dict(c)
        if c.get("check", "").startswith(("every layer density plausible", "density plausible")):
            if c.get("status") == "NOT OK":
                c["status"] = "NOTE"
            if c.get("status") == "NOTE" and "nan" in str(c.get("measured", "")):
                parts = [q.strip() for q in str(c["measured"]).split(";")]
                skipped += [q.split()[0] for q in parts if "nan" in q]
                keep = [q for q in parts if "nan" not in q]
                if keep:
                    c["measured"] = "; ".join(keep)
                else:
                    c["status"] = "ok"
        out.append(c)
    n_concern = sum(1 for c in out if c.get("status") == "NOT OK"); n_note = sum(1 for c in out if c.get("status") == "NOTE")
    if n_concern == 0:
        verdict = "PASS: all checks ok" if n_note == 0 else f"PASS: no concerns, {n_note} note(s)"
    else:
        verdict = f"REVIEW: {n_concern} concern(s)" + (f", {n_note} note(s)" if n_note else "")
    return out, verdict, (f"density check not applied to free-density layers: {', '.join(skipped)}" if skipped else None)


def draw_stack_panel(fig, out, scan, ns, title_prefix=None, extra_lines=()):
    R, main = out["result"], out["main"]
    names = list(main["names"]); v = np.asarray(main["v"], float); err = np.asarray(main["err"], float)
    p = dict(zip(names, v)); onb = {n: (n in main["at_bound"]) for n in names}; idx = {n: i for i, n in enumerate(names)}
    N = 1
    while f"t_{N + 1}" in p:
        N += 1
    lay = R["layers"]
    tt, I, lam = scan["two_theta"], scan["counts"], scan["wavelength"]
    lo, hi = R["fit_range_2theta"]; inr = (tt >= lo) & (tt <= hi)
    q = 4 * np.pi / lam * np.sin(np.deg2rad((tt - p["off"]) / 2))
    mod = ns["model_counts"](v, names, tt, lam, ns["SUBSTRATE_BETA_OVER_DELTA"], True)
    w = np.sqrt(np.maximum(I, 1)) * np.log(10)
    resid = (np.log10(np.maximum(mod, 1e-3)) - np.log10(np.maximum(I, 0.5))) * w

    # ---- right-hand column content that does not depend on the layout: summary lines (wrapped to the column width)
    lines = [f"total thickness  {R['thickness_A']:.1f} ± {R['thickness_sys_err_A']:.1f} Å (syst.; fringe estimate {R['thickness_kiessig_A']:.1f} Å)"]
    for k in range(N, 0, -1):
        L = lay[k - 1]
        if L["name"] in ns.get("_BLANK_NAMES", ()):
            lines.append(rf"{L['name']}: {L['thickness_A']:.1f} Å, SLD {p[f'sld_{k}']:.2f}$\times10^{{-5}}$ Å$^{{-2}}$ (free density)")
        else:
            lines.append(f"{L['name']}: {L['thickness_A']:.1f} Å, density {L['density_frac_of_bulk'] * 100:.0f}% of bulk")
    checks, verdict, note = fix_checks(R)
    sysfloor = rf" ($\chi^2_{{sys}}$ = {R['chi2_sys']:.2f} with a {100 * R['sys_floor']:.0f}% systematic floor)" if "chi2_sys" in R else ""
    lines.append(rf"verdict: {verdict};  reduced $\chi^2$ ={R['reduced_chi2']:.1f}{sysfloor};  residual autocorrelation (lag 1) = {R['residual_autocorr_lag1']:.2f}")
    concerns = [c for c in checks if c.get("status") == "NOT OK"]; notes_ = [c for c in checks if c.get("status") == "NOTE"]
    lines = [(ln, "0.2") for ln in lines]
    if concerns:
        lines.append(("Concerns (counted in the verdict):", "#7a2e2e")); lines += [(f"  • {c['check']}: {c['measured']}", "#7a2e2e") for c in concerns]
    if notes_:
        lines.append(("Notes (information, not counted):", "0.4")); lines += [(f"  • {c['check']}: {c['measured']}", "0.4") for c in notes_]
    if note:
        lines.append((note, "0.4"))
    lines += [(e, "0.2") for e in extra_lines]
    wrapped = [(w_, col) for ln, col in lines for w_ in (textwrap.wrap(ln, 88, subsequent_indent="   ") or [""])]
    order = " | ".join(["air"] + [lay[k - 1]["name"] for k in range(N, 0, -1)] + [R["stack"].split(" on ")[-1] + " substrate"])
    order_lines = textwrap.wrap("layers: " + order, 86, subsequent_indent="   ")
    compact = N >= 4          # many layers: one table row per layer (four parameter rows per layer would not fit); the viewer layout is kept for 1-3 layers
    eqh = (0.050, 0.050, 0.056, 0.078, 0.062) if compact else (0.060, 0.060, 0.066, 0.092, 0.075)
    nrows = 4 * N + 7; dy_full = min(0.0345, 0.46 / nrows)
    need = 0.045 + sum(eqh) + 0.032 * (len(order_lines) + 1) + 0.013 + 0.04 + ((0.030 * (N + 2) + 0.028 * 2 + 0.026 + 0.01) if compact else dy_full * nrows) + 0.012 + 0.026 * len(wrapped) + 0.02
    sc = min(1.0, 0.97 / need)          # design units assume a 13.5 x 9.5 inch figure: if the content needs more room, make the figure taller instead of cramming
    fig.set_size_inches(13.5, 9.5 / sc)

    gs = fig.add_gridspec(3, 2, width_ratios=[3, 2.2], height_ratios=[3, 3, 1.2])
    a0, a1, a2 = (fig.add_subplot(gs[i, 0]) for i in range(3))
    a0.semilogy(q, np.clip(I, 0.5, None), ".", ms=2, color="0.55", label="data")
    a0.semilogy(q[inr], mod[inr], color="crimson", lw=1.4, label="Parratt fit (layer stack)")
    a0.set_ylabel("counts"); a0.legend(fontsize=8, loc="lower left")
    sec = a0.secondary_xaxis("top", functions=(lambda x: np.degrees(2 * np.arcsin(np.clip(x * lam / (4 * np.pi), -1, 1))), lambda x: 4 * np.pi / lam * np.sin(np.radians(x) / 2)))
    sec.set_xlabel(r"$2\theta/\omega$ (deg)")
    q4 = q ** 4
    a1.semilogy(q[inr], (I * q4)[inr], ".", ms=2, color="0.55"); a1.semilogy(q[inr], (mod * q4)[inr], color="crimson", lw=1.4); a1.set_ylabel(r"$I\,q_z^4$")
    a2.plot(q[inr], resid[inr], ".", ms=2, color="0.3"); a2.axhline(0, color="crimson"); a2.set_ylabel("norm. resid."); a2.set_xlabel(r"$q_z$ ($\mathrm{\AA^{-1}}$)")

    ins = a0.inset_axes([0.52, 0.52, 0.45, 0.42])
    sld = [0.0] + [p[f"sld_{k}"] for k in range(N, 0, -1)] + [p["sld_sub"]]
    sg = [p[f"sig_{k}"] for k in range(N, 0, -1)] + [p["sig_sub"]]
    edges = np.cumsum([0.0] + [p[f"t_{k}"] for k in range(N, 0, -1)])
    z = np.linspace(-20, edges[-1] + 30, 1500); prof = np.full_like(z, sld[0])
    for j, zi in enumerate(edges):
        prof += (sld[j + 1] - sld[j]) * 0.5 * (1 + erf((z - zi) / (np.sqrt(2) * max(sg[j], 0.3))))
    ins.plot(z, prof, color="tab:blue"); ins.set_xlabel("depth (Å)", fontsize=7); ins.set_ylabel(r"SLD ($10^{-5}$ Å$^{-2}$)", fontsize=7); ins.tick_params(labelsize=7)
    for zi in edges[1:]:
        ins.axvline(zi, color="0.7", lw=0.6)

    at = fig.add_subplot(gs[:, 1]); at.axis("off"); at.set_xlim(0, 1); at.set_ylim(0, 1)
    at.text(0, 0.995, "Fitting model (Parratt recursion, Névot-Croce roughness)", fontsize=10, fontweight="bold", va="top")
    y = 0.955
    for eq, h in zip(EQ, eqh):
        at.text(0.01, y, eq, fontsize=10, va="top"); y -= h * sc
    for ol in order_lines:
        at.text(0.01, y, ol, fontsize=8, va="top", color="0.25"); y -= 0.032 * sc
    at.text(0.01, y, r"$\sigma_k$ = roughness of the top interface of layer $k$", fontsize=8, va="top", color="0.25"); y -= 0.045 * sc
    at.text(0, y, "Fitted parameters (this film)", fontsize=10, fontweight="bold", va="top"); y -= 0.04 * sc
    if not compact:
        rows = []
        for k in range(N, 0, -1):
            nm = lay[k - 1]["name"]
            rows += [(f"t_{k}", rf"d_{{{nm}}}", r"\mathrm{\AA}"), (f"sld_{k}", rf"\mathrm{{SLD}}_{{{nm}}}", r"10^{-5}\mathrm{\AA}^{-2}"),
                     (f"beta_{k}", rf"\beta_{{{nm}}}", r"10^{-5}\mathrm{\AA}^{-2}"), (f"sig_{k}", rf"\sigma_{{{nm}}}", r"\mathrm{\AA}")]
        rows += [("sld_sub", r"\mathrm{SLD}_s", r"10^{-5}\mathrm{\AA}^{-2}"), ("sig_sub", r"\sigma_s", r"\mathrm{\AA}"), ("log_scale", r"\log I_0", ""),
                 ("bg", r"I_{bg}", r"\mathrm{counts}"), ("off", r"2\theta_0", r"\mathrm{deg}"), ("fp", r"\theta_{fp}", r"\mathrm{deg}"), ("res", r"\sigma_{res}", r"\mathrm{deg}\ 2\theta")]
        for key, sym, unit in rows:
            e_ = err[idx[key]]; e_s = f" ± {e_:.2g}" if np.isfinite(e_) and e_ > 0 else ""
            nb = onb[key] and not (key.startswith(("beta_", "sig_")) or key in ("bg", "res"))
            at.text(0.01, y, f"${sym}$", fontsize=9, va="top")
            at.text(0.20, y, f"= {p[key]:.4g}{e_s}", fontsize=9, va="top", color="crimson" if nb else "k")
            at.text(0.64, y, f"${unit}$" + ("   (at bound)" if nb else "") if unit else ("(at bound)" if nb else ""), fontsize=8.5, va="top", color="0.35")
            y -= dy_full * sc
    else:
        cx = [0.01, 0.25, 0.47, 0.69, 0.85]

        def cell(key):
            e_ = err[idx[key]]; e_s = f" ± {e_:.2g}" if np.isfinite(e_) and e_ > 0 else ""
            return f"{p[key]:.4g}{e_s}", onb[key] and not (key.startswith(("beta_", "sig_")) or key in ("bg", "res"))
        for xx, h in zip(cx, ["layer", r"$d$ (Å)", "SLD", r"$\beta$", r"$\sigma$ (Å)"]):
            at.text(xx, y, h, fontsize=8.5, fontweight="bold", va="top")
        y -= 0.030 * sc
        for k in range(N, 0, -1):
            nm = lay[k - 1]["name"]; at.text(cx[0], y, nm if len(nm) <= 14 else nm[:13] + "…", fontsize=8.5, va="top")
            for xx, key in zip(cx[1:], (f"t_{k}", f"sld_{k}", f"beta_{k}", f"sig_{k}")):
                txt, nb = cell(key); at.text(xx, y, txt, fontsize=8, va="top", color="crimson" if nb else "k")
            y -= 0.030 * sc
        at.text(cx[0], y, "substrate", fontsize=8.5, va="top"); at.text(cx[1], y, "–", fontsize=8, va="top"); at.text(cx[3], y, "–", fontsize=8, va="top")
        for xx, key in ((cx[2], "sld_sub"), (cx[4], "sig_sub")):
            txt, nb = cell(key); at.text(xx, y, txt, fontsize=8, va="top", color="crimson" if nb else "k")
        y -= 0.030 * sc

        def g(key):
            return f"{p[key]:.4g} ± {err[idx[key]]:.2g}" if np.isfinite(err[idx[key]]) and err[idx[key]] > 0 else f"{p[key]:.4g}"
        at.text(0.01, y, rf"$\log I_0$ = {g('log_scale')}    $I_{{bg}}$ = {g('bg')} counts    $2\theta_0$ = {g('off')}°", fontsize=8, va="top"); y -= 0.028 * sc
        at.text(0.01, y, rf"$\theta_{{fp}}$ = {g('fp')}°" + ("  (at bound)" if onb["fp"] else "") + rf"    $\sigma_{{res}}$ = {g('res')}°", fontsize=8, va="top", color="crimson" if onb["fp"] else "k"); y -= 0.028 * sc
        at.text(0.01, y, r"SLD and $\beta$ in $10^{-5}$ Å$^{-2}$;  red = parameter at a bound", fontsize=7.5, va="top", color="0.4"); y -= 0.026 * sc
    y -= 0.012 * sc
    for w_, col in wrapped:
        at.text(0.0, y, w_, fontsize=7.8, va="top", color=col); y -= 0.026 * sc
    fig.suptitle(f"{title_prefix or R['film_id']}: XRR Parratt fit, {R['stack']}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    return fig
