# Analysis

Error analysis of the raw plaquette measurements written by `su3`, using the
gamma method ([pyerrors](https://fjosw.github.io/pyerrors/pyerrors.html)) to
account for autocorrelations.

    make analysis-setup                        # create analysis/.venv with uv
    make analyze CSV=output/plaquette_scan.csv
    make analysis-test

Options (see `python analyze.py -h`): several CSV files are combined as
replicas of one ensemble (they may differ in length, seed, thermalization
and start, but not in lattice, epsilon or sweeps per measurement); `--history [BETA ...]` and `--autocorr [BETA ...]`
write per-beta MC-history / autocorrelation plots (all betas if no value is
given); `--discard N` drops extra thermalization; `-S` sets the gamma-method
window parameter.

Output goes to `<csv stem>_analysis/`: `results.csv` (mean, error, tau_int,
window, naive error, N/tau_int, thermalization drift z-score, acceptance) and
the plots. Warnings are printed when the chain looks unthermalized or is too
short compared to tau_int.
