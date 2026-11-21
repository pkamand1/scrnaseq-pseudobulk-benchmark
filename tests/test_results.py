"""Checks on the benchmark's conclusions (run src/benchmark.py first)."""

import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
s = pd.read_csv(os.path.join(ROOT, "results", "tables", "de_benchmark_summary.csv"))


def test_pseudobulk_calls_nothing_under_the_null():
    null = s[(s.scenario == "no true DE genes") & (s.method == "pseudobulk DESeq2")]
    assert (null.mean_genes_called == 0).all()


def test_cell_level_tests_flood_false_positives():
    null = s[(s.scenario == "no true DE genes") & s.method.str.startswith("cell-level")]
    assert (null.mean_genes_called > 100).all()


def test_pseudobulk_fdp_near_target_and_power_grows_with_donors():
    pb = s[(s.scenario == "100 true DE genes") & (s.method == "pseudobulk DESeq2")].sort_values("donors_per_group")
    assert (pb.mean_fdp < 0.10).all()
    assert pb.mean_power.is_monotonic_increasing
