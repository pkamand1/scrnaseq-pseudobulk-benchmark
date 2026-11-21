"""Cell-level versus pseudobulk differential expression on simulated
multi-donor single-cell data, plus a clustering resolution sweep.

Why this matters: cells from the same donor are not independent. Treating
each cell as a replicate (pseudoreplication) ignores donor-to-donor
variability, which inflates false discoveries; summing counts per donor and
cell type (pseudobulk) and testing across donors respects the design.

Simulation: 2,000 genes, five cell types with marker genes, donor-specific
expression shifts (biological variation between people), varying library
sizes, and negative binomial counts. In T cells, 100 genes truly change with
condition (log2 fold change of ±1); every other gene is null. Each scenario
is repeated over independent simulated studies.

Usage: python src/benchmark.py [--reps 8]
"""

import argparse
import os
import warnings

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G, N_DE, CELL_TYPES = 2000, 100, ["T cell", "B cell", "NK cell", "Monocyte", "Dendritic"]
TYPE_P = [0.40, 0.20, 0.12, 0.22, 0.06]


def simulate(n_donors_per_group, cells_per_donor, seed, de=True, donor_sd=0.25):
    rng = np.random.default_rng(seed)
    base = rng.lognormal(np.log(0.4), 1.2, G)                    # mean expression per gene
    type_lfc = np.zeros((len(CELL_TYPES), G))
    for t in range(len(CELL_TYPES)):                            # 40 markers per type, 2.6-fold up
        type_lfc[t, rng.choice(G, 40, replace=False)] = 1.4
    de_genes = rng.choice(np.arange(G), N_DE, replace=False) if de else np.array([], int)
    cond_lfc = np.zeros(G)
    cond_lfc[de_genes] = rng.choice([-1.0, 1.0], len(de_genes))
    rows, mats = [], []
    for cond in (0, 1):
        for d in range(n_donors_per_group):
            donor = f"{'ctrl' if cond == 0 else 'stim'}_{d + 1:02d}"
            donor_eff = rng.normal(0, donor_sd, G)                # person-specific log shift
            types = rng.choice(len(CELL_TYPES), cells_per_donor, p=TYPE_P)
            lib = rng.lognormal(0, 0.35, cells_per_donor)
            log_mu = (np.log(base)[None, :] + np.log(2) * type_lfc[types] + donor_eff[None, :]
                      + np.log(2) * cond * cond_lfc[None, :] * (types == 0)[:, None])   # condition acts in T cells
            mu = np.exp(log_mu) * lib[:, None]
            counts = rng.poisson(rng.gamma(5.0, mu / 5.0))        # negative binomial, dispersion 0.2
            mats.append(counts.astype(np.float32))
            rows += [{"donor": donor, "condition": "stim" if cond else "ctrl", "cell_type": CELL_TYPES[t]} for t in types]
    a = ad.AnnData(np.vstack(mats), obs=pd.DataFrame(rows))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.var_names = [f"g{i:04d}" for i in range(G)]
    truth = np.zeros(G, bool)
    truth[de_genes] = True
    return a, truth


def cell_level(a, method):
    t = a[a.obs.cell_type == "T cell"].copy()
    sc.pp.normalize_total(t, target_sum=1e4)
    sc.pp.log1p(t)
    sc.tl.rank_genes_groups(t, "condition", groups=["stim"], reference="ctrl", method=method)
    r = sc.get.rank_genes_groups_df(t, "stim").set_index("names").reindex(t.var_names)
    return r.pvals.fillna(1).values


def pseudobulk(a):
    t = a[a.obs.cell_type == "T cell"]
    donors = t.obs.donor.unique()
    counts = pd.DataFrame(np.vstack([t[t.obs.donor == d].X.sum(axis=0) for d in donors]).astype(int),
                          index=donors, columns=t.var_names)
    meta = pd.DataFrame({"condition": [d.split("_")[0] for d in donors]}, index=donors)
    keep = counts.sum(axis=0) >= 10
    dds = DeseqDataSet(counts=counts.loc[:, keep], metadata=meta, design="~condition", quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "stim", "ctrl"], quiet=True)
    st.summary()
    p = pd.Series(1.0, index=t.var_names)
    p[keep[keep].index] = st.results_df.pvalue.fillna(1).values
    return p.values


def evaluate(p, truth, alpha=0.05):
    rej = multipletests(p, alpha=alpha, method="fdr_bh")[0]
    fdp = (rej & ~truth).sum() / max(rej.sum(), 1)
    power = (rej & truth).sum() / truth.sum() if truth.any() else np.nan
    return rej.sum(), fdp, power, (rej & ~truth).sum()


def clustering_sweep(seed=1):
    a, _ = simulate(6, 300, seed)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=500)
    sc.pp.pca(a, n_comps=30, mask_var="highly_variable")
    sc.pp.neighbors(a, n_neighbors=15)
    rows = []
    for res in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]:
        sc.tl.leiden(a, resolution=res, key_added="l", flavor="igraph", n_iterations=2, random_state=0)
        rows.append({"resolution": res, "clusters": a.obs.l.nunique(),
                     "adjusted_rand_index": round(adjusted_rand_score(a.obs.cell_type, a.obs.l), 3)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=8)
    args = ap.parse_args()
    rows = []
    for n_donors in (3, 6, 10):
        for rep in range(args.reps):
            for de in (True, False):
                a, truth = simulate(n_donors, 300, seed=1000 * n_donors + 10 * rep + int(de), de=de)
                for name, f in [("cell-level Wilcoxon", lambda x: cell_level(x, "wilcoxon")),
                                ("cell-level t-test", lambda x: cell_level(x, "t-test")),
                                ("pseudobulk DESeq2", pseudobulk)]:
                    n_rej, fdp, power, n_false = evaluate(f(a), truth)
                    rows.append({"donors_per_group": n_donors, "rep": rep,
                                 "scenario": "100 true DE genes" if de else "no true DE genes",
                                 "method": name, "genes_called": n_rej, "false_discoveries": n_false,
                                 "false_discovery_proportion": fdp, "power": power})
            print(f"donors {n_donors}, replicate {rep + 1} done", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(ROOT, "results", "tables", "de_benchmark_all_runs.csv"), index=False)
    summ = (res.groupby(["scenario", "donors_per_group", "method"])
            .agg(mean_genes_called=("genes_called", "mean"), mean_false_discoveries=("false_discoveries", "mean"),
                 mean_fdp=("false_discovery_proportion", "mean"), mean_power=("power", "mean"))
            .round(3).reset_index())
    summ.to_csv(os.path.join(ROOT, "results", "tables", "de_benchmark_summary.csv"), index=False)
    clus = clustering_sweep()
    clus.to_csv(os.path.join(ROOT, "results", "tables", "clustering_resolution.csv"), index=False)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    colors = {"cell-level Wilcoxon": "#d62728", "cell-level t-test": "#ff7f0e", "pseudobulk DESeq2": "#1f77b4"}
    s1 = summ[summ.scenario == "100 true DE genes"]
    s0 = summ[summ.scenario == "no true DE genes"]
    for m, g in s1.groupby("method"):
        ax[0].plot(g.donors_per_group, g.mean_fdp, "o-", color=colors[m], label=m)
        ax[1].plot(g.donors_per_group, g.mean_power, "o-", color=colors[m], label=m)
    for m, g in s0.groupby("method"):
        ax[2].plot(g.donors_per_group, g.mean_false_discoveries, "o-", color=colors[m], label=m)
    ax[0].axhline(0.05, ls=":", color="k")
    ax[0].set_title("False discovery proportion (target 5%)")
    ax[1].set_title("Power: share of 100 true DE genes found")
    ax[2].set_title("Genes called when nothing changes")
    ax[2].set_yscale("symlog")
    for x in ax:
        x.set_xlabel("Donors per condition")
        x.set_xticks([3, 6, 10])
    ax[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "results", "figures", "pseudobulk_vs_cell_level.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(clus.resolution, clus.adjusted_rand_index, "o-")
    for _, r in clus.iterrows():
        ax.annotate(f"{int(r.clusters)}", (r.resolution, r.adjusted_rand_index), xytext=(4, 4),
                    textcoords="offset points", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("Leiden resolution (labels: number of clusters)")
    ax.set_ylabel("Adjusted Rand index vs true cell types")
    ax.set_title("Clustering resolution sweep (5 true cell types)")
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "results", "figures", "clustering_resolution.png"), dpi=150)
    plt.close(fig)
    print(summ.to_string(index=False))
    print(clus.to_string(index=False))


if __name__ == "__main__":
    main()
