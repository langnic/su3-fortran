#!/usr/bin/env python
"""Error analysis of the raw plaquette measurements written by the su3 program.

Uses the gamma method (Wolff 2004; pyerrors) to estimate statistical errors
including autocorrelations, and integrated autocorrelation times tau_int.

    python analyze.py ../output/plaquette_scan.csv
    python analyze.py run_a.csv run_b.csv        # independent runs = replicas of one ensemble
    python analyze.py out.csv --history          # MC history plot for every beta
    python analyze.py out.csv --history 0.5 6.0 --autocorr 6.0

Output (default: <csv stem>_analysis/ next to the first input file):
    results.csv              one row per beta: mean, error, tau_int, diagnostics
    plaquette_vs_beta.*      <P> with gamma-method errors
    tauint_vs_beta.*         integrated autocorrelation time
    acceptance_vs_beta.*     Metropolis acceptance (useful for tuning epsilon)
    history_beta_*.*         only with --history
    autocorr_beta_*.*        only with --autocorr
"""
import argparse
import csv
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyerrors as pe  # noqa: E402

BETA_DECIMALS = 4  # betas are written with 4 decimals

# Header keys that may differ between replicas without making them incompatible: run length,
# start and thermalization, bookkeeping. Everything else (lattice, epsilon, sweeps per
# measurement, ...) changes the Markov chain itself, and the gamma method pools the
# autocorrelation of all replicas, so those must match.
REPLICA_KEYS = {"seed", "config", "images", "measurements", "thermalization_sweeps",
                "hot_start", "reunitarize_interval"}


# --------------------------------------------------------------------------- input

def read_measurements(path):
    """Read one raw-data file. Returns (meta dict, {beta: dict of column arrays})."""
    meta = {}
    columns = None
    rows = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "=" in line:
                    key, _, value = line[1:].partition("=")
                    meta[key.strip()] = value.strip()
                continue
            fields = [x.strip() for x in line.split(",")]
            if columns is None:
                columns = fields
                missing = {"beta", "plaquette"} - set(columns)
                if missing:
                    raise ValueError(f"{path}: missing column(s) {sorted(missing)}")
                continue
            if len(fields) != len(columns):
                raise ValueError(f"{path}:{lineno}: expected {len(columns)} fields, got {len(fields)}")
            rows.append([float(x) for x in fields])
    if columns is None or not rows:
        raise ValueError(f"{path}: no data")

    data = np.array(rows)
    by_beta = {}
    beta_col = data[:, columns.index("beta")]
    for beta in dict.fromkeys(np.round(beta_col, BETA_DECIMALS)):   # keeps file order
        sel = np.isclose(beta_col, beta, atol=0.5 * 10.0**-BETA_DECIMALS)
        by_beta[float(beta)] = {name: data[sel, i] for i, name in enumerate(columns) if name != "beta"}
    return meta, by_beta


def check_compatible(metas, paths):
    """Replicas must come from the same physical setup."""
    ref = metas[0]
    for meta, path in zip(metas[1:], paths[1:]):
        for key in sorted((set(ref) | set(meta)) - REPLICA_KEYS):
            if ref.get(key) != meta.get(key):
                raise ValueError(
                    f"{path} is not compatible with {paths[0]}: '{key}' differs "
                    f"({ref.get(key)!r} vs {meta.get(key)!r})")


# --------------------------------------------------------------------------- analysis

def ensemble_name(beta):
    return f"beta_{beta:.{BETA_DECIMALS}f}"


def make_obs(replicas, name, discard=0):
    """pyerrors observable from a list of 1d sample arrays (one per replica)."""
    samples = [np.asarray(r)[discard:] for r in replicas]
    if any(len(s) < 5 for s in samples):
        raise ValueError(f"{name}: need at least 5 measurements per replica after discarding {discard}")
    names = [name] if len(samples) == 1 else [f"{name}|r{i}" for i in range(len(samples))]
    return pe.Obs(samples, names)


def _drift_single(r, name, S):
    """z-score of (first half - second half) for one chain; nan if too short to split."""
    h = len(r) // 2
    if h < 5:
        return np.nan
    a = pe.Obs([r[:h]], [f"{name}_a"])
    b = pe.Obs([r[h:2 * h]], [f"{name}_b"])
    a.gamma_method(S=S)
    b.gamma_method(S=S)
    err = np.hypot(a.dvalue, b.dvalue)
    diff = a.value - b.value
    if err == 0:
        return 0.0 if diff == 0 else np.inf
    return diff / err


def drift_z(replicas, name, discard, S):
    """Compare first and second half of each chain: z-score of the difference.
    Large |z| means a chain has not equilibrated (or drifts). Computed per replica, so that
    one badly thermalized run is not hidden by merging it with good ones; the value with the
    largest |z| is returned. Replicas too short to split are skipped (nan if none is left)."""
    zs = [_drift_single(np.asarray(r)[discard:], f"{name}_r{i}", S) for i, r in enumerate(replicas)]
    zs = [z for z in zs if not np.isnan(z)]
    return max(zs, key=abs) if zs else np.nan


def analyze_beta(beta, replicas, acceptance, sweeps_per_meas, S=2.0, discard=0):
    """Gamma-method analysis of the plaquette at one beta. Returns (result dict, Obs)."""
    name = ensemble_name(beta)
    obs = make_obs(replicas, name, discard)
    obs.gamma_method(S=S)

    n = int(sum(len(r) - discard for r in replicas))
    tau = float(obs.e_tauint[name])
    flat = np.concatenate([np.asarray(r)[discard:] for r in replicas])
    result = {
        "beta": beta,
        "n_meas": n,
        "plaquette": float(obs.value),
        "error": float(obs.dvalue),
        "tau_int": tau,                                    # in units of measurements
        "dtau_int": float(obs.e_dtauint[name]),
        "tau_int_sweeps": tau * sweeps_per_meas,
        "window": int(obs.e_windowsize[name]),
        "naive_error": float(np.std(flat, ddof=1) / np.sqrt(n)),
        "n_over_tau": n / tau,
        "drift_z": float(drift_z(replicas, name, discard, S)),
        "acceptance": float(np.mean(np.concatenate([np.asarray(a)[discard:] for a in acceptance]))),
    }
    return result, obs


def warnings_for(r):
    """Human-readable warnings for one result row."""
    out = []
    if abs(r["drift_z"]) > 3:
        out.append(f"beta={r['beta']:g}: first and second half of the chain differ by {r['drift_z']:.1f} sigma "
                   "-> not thermalized? (try --discard or more thermalization sweeps)")
    if r["n_over_tau"] < 50:
        out.append(f"beta={r['beta']:g}: only {r['n_over_tau']:.0f} measurements per tau_int "
                   "-> error estimate unreliable, run longer")
    return out


# --------------------------------------------------------------------------- output

RESULT_FIELDS = ["beta", "n_meas", "plaquette", "error", "tau_int", "dtau_int", "tau_int_sweeps",
                 "window", "naive_error", "n_over_tau", "drift_z", "acceptance"]


def write_results(results, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(RESULT_FIELDS)
        for r in results:
            w.writerow([f"{r[k]:.{BETA_DECIMALS}f}" if k == "beta" else
                        (r[k] if k in ("n_meas", "window") else f"{r[k]:.8g}") for k in RESULT_FIELDS])


def lattice_label(meta):
    dims = meta.get("lattice (nx ny nz nt)")
    if not dims:
        return ""
    d = dims.split()
    return f"${d[0]}\\times{d[1]}\\times{d[2]}\\times{d[3]}$ lattice" if len(d) == 4 else ""


def save(fig, base, formats):
    for fmt in formats:
        fig.savefig(f"{base}.{fmt}", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_plaquette(results, meta, base, formats, reference=True):
    beta = np.array([r["beta"] for r in results])
    fig, ax = plt.subplots(figsize=(6, 4))
    if reference and len(beta) > 0:
        x = np.linspace(0, max(beta.max(), 1e-9), 100)
        ax.plot(x, x / 18, color="0.6", ls="--", lw=1, label=r"strong coupling, $\beta/18$")
    ax.errorbar(beta, [r["plaquette"] for r in results], yerr=[r["error"] for r in results],
                fmt="o", ms=3, capsize=2, label=lattice_label(meta) or None)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\langle P\rangle$")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    save(fig, base, formats)


def plot_tauint(results, meta, base, formats):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar([r["beta"] for r in results], [r["tau_int_sweeps"] for r in results],
                yerr=[r["dtau_int"] * float(meta.get("sweeps_per_measurement", 1)) for r in results],
                fmt="o", ms=4, capsize=2)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\tau_\mathrm{int}$ [sweeps]")
    ax.grid(alpha=0.3)
    save(fig, base, formats)


def plot_acceptance(results, base, formats):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot([r["beta"] for r in results], [r["acceptance"] for r in results], "o-", ms=4)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel("Metropolis acceptance")
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    save(fig, base, formats)


def plot_history(beta, replicas, obs, base, formats, discard=0):
    """MC history of the plaquette with mean and error band, plus histogram."""
    fig, (ax, axh) = plt.subplots(1, 2, figsize=(8, 3.5), sharey=True,
                                  gridspec_kw={"width_ratios": [4, 1], "wspace": 0.05})
    for i, rep in enumerate(replicas):
        y = np.asarray(rep)[discard:]
        ax.plot(np.arange(discard, discard + len(y)), y, ".", ms=3, label=f"replica {i}" if len(replicas) > 1 else None)
    ax.axhline(obs.value, color="k", lw=1)
    ax.axhspan(obs.value - obs.dvalue, obs.value + obs.dvalue, color="k", alpha=0.15)
    ax.set_xlabel("measurement")
    ax.set_ylabel(r"$P$")
    ax.set_title(rf"$\beta = {beta:g}$:  $\langle P\rangle = {obs.value:.5f} \pm {obs.dvalue:.5f}$", fontsize=10)
    if len(replicas) > 1:
        ax.legend(frameon=False, fontsize=8)
    axh.hist(np.concatenate([np.asarray(r)[discard:] for r in replicas]), bins=20, orientation="horizontal", alpha=0.7)
    axh.set_xticks([])
    save(fig, base, formats)


def plot_autocorr(beta, obs, base, formats):
    """Normalized autocorrelation function and tau_int(W) with the chosen summation window."""
    name = ensemble_name(beta)
    rho, drho = np.asarray(obs.e_rho[name]), np.asarray(obs.e_drho[name])
    tau_w, dtau_w = np.asarray(obs.e_n_tauint[name]), np.asarray(obs.e_n_dtauint[name])
    w = int(obs.e_windowsize[name])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.5))
    # the tail of both functions is pure noise: show only a few window lengths
    tmax = min(len(rho), max(3 * w, 30))
    t = np.arange(tmax)
    a1.plot(t, rho[:tmax], "-", lw=1)
    a1.fill_between(t, (rho - drho)[:tmax], (rho + drho)[:tmax], alpha=0.25)
    a1.axhline(0, color="k", lw=0.5)
    a1.axvline(w, color="C3", ls="--", lw=1)
    a1.set_xlabel("$t$ [measurements]")
    a1.set_ylabel(r"$\rho(t)$")
    tmax_w = min(len(tau_w), tmax)
    tw = np.arange(tmax_w)
    a2.plot(tw, tau_w[:tmax_w], "-", lw=1)
    a2.fill_between(tw, (tau_w - dtau_w)[:tmax_w], (tau_w + dtau_w)[:tmax_w], alpha=0.25)
    a2.axvline(w, color="C3", ls="--", lw=1, label=f"window $W={w}$")
    a2.set_xlabel("$W$ [measurements]")
    a2.set_ylabel(r"$\tau_\mathrm{int}(W)$")
    a2.legend(frameon=False)
    fig.suptitle(rf"$\beta = {beta:g}$", fontsize=10)
    fig.tight_layout()
    save(fig, base, formats)


# --------------------------------------------------------------------------- driver

def select_betas(requested, available):
    """Map requested betas onto the available grid. [] (flag without values) means all."""
    if requested is None:
        return []
    if len(requested) == 0:
        return list(available)
    chosen = []
    for b in requested:
        nearest = min(available, key=lambda a: abs(a - b))
        if abs(nearest - b) > 10.0**-BETA_DECIMALS * 10:
            raise ValueError(f"beta = {b} not found in data (available: {', '.join(f'{a:g}' for a in available)})")
        chosen.append(nearest)
    return chosen


def run(paths, out_dir=None, S=2.0, discard=0, history=None, autocorr=None, formats=("png",), reference=True):
    paths = [Path(p) for p in paths]
    loaded = [read_measurements(p) for p in paths]
    metas = [m for m, _ in loaded]
    check_compatible(metas, [str(p) for p in paths])
    meta = metas[0]
    sweeps_per_meas = int(meta.get("sweeps_per_measurement", 1))

    betas = sorted(set().union(*[set(d) for _, d in loaded]))
    for (_, d), p in zip(loaded, paths):
        if set(d) != set(betas):
            raise ValueError(f"{p} has a different set of beta values than the other files")

    if out_dir is None:
        out_dir = paths[0].parent / f"{paths[0].stem}_analysis"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if len(paths) > 1:
        lengths = ", ".join(f"{p.name}: {len(next(iter(d.values()))['plaquette'])}" for (_, d), p in zip(loaded, paths))
        print(f"merging {len(paths)} replicas, measurements per beta ({lengths})")

    want_history = select_betas(history, betas)
    want_autocorr = select_betas(autocorr, betas)

    results, all_warnings = [], []
    for beta in betas:
        replicas = [d[beta]["plaquette"] for _, d in loaded]
        accept = [d[beta].get("acceptance", np.array([np.nan])) for _, d in loaded]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)   # e.g. zero variance (cold start)
            r, obs = analyze_beta(beta, replicas, accept, sweeps_per_meas, S=S, discard=discard)
        results.append(r)
        all_warnings += warnings_for(r)
        if beta in want_history:
            plot_history(beta, replicas, obs, out_dir / f"history_beta_{beta:.{BETA_DECIMALS}f}", formats, discard)
        if beta in want_autocorr:
            plot_autocorr(beta, obs, out_dir / f"autocorr_beta_{beta:.{BETA_DECIMALS}f}", formats)

    write_results(results, out_dir / "results.csv")
    plot_plaquette(results, meta, out_dir / "plaquette_vs_beta", formats, reference)
    plot_tauint(results, meta, out_dir / "tauint_vs_beta", formats)
    plot_acceptance(results, out_dir / "acceptance_vs_beta", formats)
    return results, all_warnings, out_dir


def print_table(results):
    print(f"{'beta':>7} {'N':>6} {'<P>':>10} {'error':>10} {'naive':>10} {'tau_int':>8} {'+-':>6} {'W':>4} {'N/tau':>7} {'drift z':>8}")
    for r in results:
        print(f"{r['beta']:7.3f} {r['n_meas']:6d} {r['plaquette']:10.6f} {r['error']:10.2e} {r['naive_error']:10.2e} "
              f"{r['tau_int']:8.2f} {r['dtau_int']:6.2f} {r['window']:4d} {r['n_over_tau']:7.0f} {r['drift_z']:8.2f}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+", help="raw measurement file(s); several files are treated as replicas "
                    "(they may differ in length, seed, thermalization and start, but not in lattice or algorithm settings)")
    ap.add_argument("-o", "--out", help="output directory (default: <csv stem>_analysis next to the input)")
    ap.add_argument("-S", type=float, default=2.0, help="gamma-method window parameter (default 2.0)")
    ap.add_argument("--discard", type=int, default=0, metavar="N",
                    help="discard the first N measurements of each chain (default 0)")
    ap.add_argument("--history", nargs="*", type=float, metavar="BETA",
                    help="write MC history plots; without values for every beta")
    ap.add_argument("--autocorr", nargs="*", type=float, metavar="BETA",
                    help="write autocorrelation function / tau_int(W) plots; without values for every beta")
    ap.add_argument("--format", nargs="+", default=["png"], choices=["png", "pdf", "svg"],
                    help="figure file format(s) (default png)")
    ap.add_argument("--no-reference", action="store_true", help="omit the strong-coupling line from the plaquette plot")
    args = ap.parse_args(argv)

    try:
        results, warns, out_dir = run(args.csv, args.out, args.S, args.discard, args.history, args.autocorr,
                                      tuple(args.format), not args.no_reference)
    except (ValueError, OSError) as e:
        sys.exit(f"error: {e}")

    print_table(results)
    for w in warns:
        print("WARNING:", w, file=sys.stderr)
    print(f"\nresults and plots written to {out_dir}/")


if __name__ == "__main__":
    main()
