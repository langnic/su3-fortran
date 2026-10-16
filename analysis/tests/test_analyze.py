"""Tests for the analysis script, using synthetic chains with known statistics."""
import numpy as np
import pytest

import analyze


def ar1(n, a, rng, mean=0.0, sigma_stat=1.0):
    """AR(1) chain: x_i = a x_{i-1} + noise. Exact tau_int = (1+a) / (2 (1-a))."""
    x = np.empty(n)
    x[0] = rng.normal()
    for i in range(1, n):
        x[i] = a * x[i - 1] + rng.normal()
    return mean + sigma_stat * x / x.std()


def tau_exact(a):
    return (1 + a) / (2 * (1 - a))


def write_csv(path, chains, header=None, acceptance=0.5):
    """chains: {beta: samples}"""
    header = {"seed": "1", "lattice (nx ny nz nt)": "4 4 4 4", "sweeps_per_measurement": "2", **(header or {})}
    with open(path, "w") as f:
        f.write("# su3-fortran raw plaquette measurements\n")
        for k, v in header.items():
            f.write(f"# {k} = {v}\n")
        f.write("beta,measurement,plaquette,acceptance\n")
        for beta, x in chains.items():
            for i, p in enumerate(x, 1):
                f.write(f"{beta:8.4f},{i:6d},{p:16.8e},{acceptance:8.5f}\n")


def test_recovers_tau_int_and_error():
    rng = np.random.default_rng(10)
    a, n = 0.8, 20000
    x = ar1(n, a, rng, mean=0.3, sigma_stat=0.01)
    res, _ = analyze.analyze_beta(1.0, [x], [np.full(n, 0.5)], sweeps_per_meas=2)
    assert res["tau_int"] == pytest.approx(tau_exact(a), rel=0.15)
    assert res["tau_int_sweeps"] == pytest.approx(2 * res["tau_int"])
    exact_err = 0.01 / np.sqrt(n) * np.sqrt(2 * tau_exact(a))
    assert res["error"] == pytest.approx(exact_err, rel=0.15)
    # the whole point: the naive error is too small by sqrt(2 tau_int)
    assert res["error"] / res["naive_error"] == pytest.approx(np.sqrt(2 * res["tau_int"]), rel=0.05)


def test_uncorrelated_data_has_tau_one_half():
    rng = np.random.default_rng(11)
    x = rng.normal(size=20000)
    res, _ = analyze.analyze_beta(1.0, [x], [np.zeros(len(x))], 1)
    assert res["tau_int"] == pytest.approx(0.5, abs=0.1)
    assert res["error"] == pytest.approx(res["naive_error"], rel=0.1)


def test_constant_chain_does_not_crash():
    res, _ = analyze.analyze_beta(0.5, [np.ones(100)], [np.ones(100)], 1)
    assert res["error"] == 0.0 and res["plaquette"] == 1.0
    assert res["drift_z"] == 0.0


def test_drift_detected_for_unthermalized_chain():
    rng = np.random.default_rng(12)
    x = ar1(4000, 0.5, rng, sigma_stat=0.01) + np.linspace(0.1, 0, 4000)   # relaxation towards 0
    z = analyze.drift_z([x], "t", 0, 2.0)
    assert abs(z) > 3
    # discarding the transient reduces it
    assert abs(analyze.drift_z([x], "t", 2000, 2.0)) < abs(z)
    # a stationary chain is fine
    y = ar1(4000, 0.5, rng)
    assert abs(analyze.drift_z([y], "t", 0, 2.0)) < 3


def test_too_few_measurements_is_a_clear_error():
    with pytest.raises(ValueError, match="at least 5"):
        analyze.make_obs([np.arange(3.0)], "x")


def test_read_measurements_header_and_grouping(tmp_path):
    rng = np.random.default_rng(13)
    p = tmp_path / "a.csv"
    write_csv(p, {0.5: rng.normal(size=30), 1.5: rng.normal(size=30)}, header={"epsilon": ".30000"})
    meta, data = analyze.read_measurements(p)
    assert meta["epsilon"] == ".30000" and meta["sweeps_per_measurement"] == "2"
    assert list(data) == [0.5, 1.5]
    assert len(data[0.5]["plaquette"]) == 30 and data[0.5]["acceptance"][0] == 0.5


def test_read_measurements_rejects_bad_files(tmp_path):
    (tmp_path / "empty.csv").write_text("# only header\n")
    with pytest.raises(ValueError, match="no data"):
        analyze.read_measurements(tmp_path / "empty.csv")
    (tmp_path / "ragged.csv").write_text("beta,measurement,plaquette\n1.0,1\n")
    with pytest.raises(ValueError, match="expected 3 fields"):
        analyze.read_measurements(tmp_path / "ragged.csv")
    (tmp_path / "nocol.csv").write_text("beta,measurement\n1.0,1\n")
    with pytest.raises(ValueError, match="plaquette"):
        analyze.read_measurements(tmp_path / "nocol.csv")


def test_replicas_are_combined_and_incompatible_ones_rejected(tmp_path):
    rng = np.random.default_rng(14)
    a, b, c = tmp_path / "a.csv", tmp_path / "b.csv", tmp_path / "c.csv"
    write_csv(a, {1.0: ar1(2000, 0.5, rng, 0.2, 0.01)}, header={"seed": "1"})
    write_csv(b, {1.0: ar1(2000, 0.5, rng, 0.2, 0.01)}, header={"seed": "2"})
    write_csv(c, {1.0: ar1(2000, 0.5, rng, 0.2, 0.01)}, header={"seed": "3", "epsilon": "0.9"})
    results, _, _ = analyze.run([a, b], out_dir=tmp_path / "out")
    assert results[0]["n_meas"] == 4000
    assert results[0]["plaquette"] == pytest.approx(0.2, abs=0.002)
    with pytest.raises(ValueError, match="epsilon"):
        analyze.run([a, c], out_dir=tmp_path / "out2")


def test_end_to_end_outputs_and_optional_plots(tmp_path):
    rng = np.random.default_rng(15)
    chains = {b: ar1(1000, 0.6, rng, b / 18, 0.005) for b in (0.5, 1.0, 1.5)}
    p = tmp_path / "scan.csv"
    write_csv(p, chains)

    # default: summary plots only
    results, warns, out = analyze.run([p])
    assert out == tmp_path / "scan_analysis"
    for name in ("results.csv", "plaquette_vs_beta.png", "tauint_vs_beta.png", "acceptance_vs_beta.png"):
        assert (out / name).stat().st_size > 0
    assert not list(out.glob("history_*")) and not list(out.glob("autocorr_*"))
    assert [r["beta"] for r in results] == [0.5, 1.0, 1.5]
    lines = (out / "results.csv").read_text().splitlines()
    assert lines[0].startswith("beta,n_meas,plaquette,error") and len(lines) == 4

    # requested: history for one beta, autocorr for all, pdf as well
    _, _, out2 = analyze.run([p], out_dir=tmp_path / "o2", history=[1.0], autocorr=[], formats=("png", "pdf"))
    assert [f.name for f in sorted(out2.glob("history_*"))] == ["history_beta_1.0000.pdf", "history_beta_1.0000.png"]
    assert len(list(out2.glob("autocorr_*.png"))) == 3


def test_unknown_beta_request_is_an_error(tmp_path):
    rng = np.random.default_rng(16)
    p = tmp_path / "scan.csv"
    write_csv(p, {1.0: rng.normal(size=100)})
    with pytest.raises(ValueError, match="not found"):
        analyze.run([p], out_dir=tmp_path / "o", history=[2.0])


def test_cli_exit_on_missing_file(tmp_path):
    with pytest.raises(SystemExit):
        analyze.main([str(tmp_path / "nope.csv")])


def test_replicas_of_different_length_are_merged(tmp_path):
    rng = np.random.default_rng(20)
    a, b = 0.7, 0.7
    short = ar1(100, a, rng, 0.2, 0.01)
    long_ = ar1(300, b, rng, 0.2, 0.01)
    f1, f2 = tmp_path / "short.csv", tmp_path / "long.csv"
    # runs of different length, thermalization and start are fine
    write_csv(f1, {1.0: short}, header={"seed": "1", "measurements": "100", "thermalization_sweeps": "50", "hot_start": "F"})
    write_csv(f2, {1.0: long_}, header={"seed": "2", "measurements": "300", "thermalization_sweeps": "200", "hot_start": "T"})
    results, _, _ = analyze.run([f1, f2], out_dir=tmp_path / "out", history=[1.0], autocorr=[1.0])
    r = results[0]
    assert r["n_meas"] == 400
    # pyerrors pools the replicas: the central value is the mean of all samples
    assert r["plaquette"] == pytest.approx(np.concatenate([short, long_]).mean(), abs=1e-9)   # the csv stores 9 significant digits
    assert r["error"] > r["naive_error"]
    assert (tmp_path / "out" / "history_beta_1.0000.png").exists()
    assert (tmp_path / "out" / "autocorr_beta_1.0000.png").exists()


def test_replicas_with_different_chain_settings_are_rejected(tmp_path):
    rng = np.random.default_rng(21)
    base = tmp_path / "base.csv"
    write_csv(base, {1.0: rng.normal(size=100)})
    for key, value in [("sweeps_per_measurement", "5"), ("epsilon", "0.9"), ("lattice (nx ny nz nt)", "8 8 8 8")]:
        other = tmp_path / "other.csv"
        write_csv(other, {1.0: rng.normal(size=100)}, header={key: value})
        with pytest.raises(ValueError, match="differs"):
            analyze.run([base, other], out_dir=tmp_path / "o")


def test_drift_ignores_replicas_too_short_to_split():
    rng = np.random.default_rng(22)
    z = analyze.drift_z([rng.normal(size=500), rng.normal(size=6)], "t", 0, 2.0)
    assert np.isfinite(z)
    assert np.isnan(analyze.drift_z([rng.normal(size=6)], "t", 0, 2.0))


def test_drift_is_per_replica_so_a_bad_run_is_not_diluted():
    rng = np.random.default_rng(23)
    good = ar1(2000, 0.5, rng, 0.0, 0.01)
    bad = ar1(300, 0.5, rng, 0.0, 0.01) + np.linspace(0.05, 0, 300)   # still relaxing
    z_bad = analyze.drift_z([bad], "t", 0, 2.0)
    assert abs(z_bad) > 3
    assert analyze.drift_z([good, bad], "t", 0, 2.0) == pytest.approx(z_bad)
    assert abs(analyze.drift_z([good], "t", 0, 2.0)) < 3
