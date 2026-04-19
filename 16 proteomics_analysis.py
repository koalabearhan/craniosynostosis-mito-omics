
from pathlib import Path
import math
import re
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact, binomtest, ttest_ind, norm
from statsmodels.stats.multitest import multipletests

BASE = Path(__file__).resolve().parents[1] / "input_files"

def load_inputs(base: Path):
    vol = pd.read_excel(base / "1.data.xlsx")
    prot = pd.read_excel(base / "PROTEOMICS result.xlsx", sheet_name="Protein quantitation")
    mito = pd.read_excel(base / "Human Mitocarta 3.0.xlsx")
    mitopath = pd.read_excel(base / "mitopathways.xlsx")
    trans = pd.read_csv(base / "differential_expression_results.csv")
    # first column is gene symbol
    trans = trans.rename(columns={trans.columns[0]: "GeneSymbol"})
    return vol, prot, mito, mitopath, trans

def first_symbol(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return np.nan
    # split vendor multi-symbol fields like "GeneA;GeneB"
    s = re.split(r"[;,/]", s)[0].strip()
    return s if s else np.nan

def parse_gene_list(x):
    if pd.isna(x):
        return []
    return [g.strip().upper() for g in str(x).split(",") if g and str(g).lower() != "nan"]

def main():
    base = BASE
    vol, prot, mito, mitopath, trans = load_inputs(base)

    # Normalize gene symbol columns
    prot["GeneSymbol"] = prot["Gene Name"].apply(first_symbol)
    vol = vol.merge(prot[["Accession", "GeneSymbol"]], left_on="Protein", right_on="Accession", how="left")
    vol["GeneUpper"] = vol["GeneSymbol"].astype(str).str.upper()
    vol = vol[vol["GeneUpper"] != "NAN"].copy()

    trans["GeneUpper"] = trans["GeneSymbol"].astype(str).str.upper()

    mito["GeneUpper"] = mito["Symbol"].astype(str).str.upper()
    mito_symbols = set(mito["GeneUpper"].dropna())
    mitopath = mitopath.rename(columns={2: "ID"})
    mitopath["GeneSet"] = mitopath["Genes"].apply(parse_gene_list)

    # Current candidate families used in the revised manuscript
    mitodeg24 = ['SLC25A21','CYP27A1','GSR','MRPS30','ACAA2','CA5A','LIG3','ALDH2','GPT2','PNPT1',
                 'ME3','MTHFD2','GLS2','PINK1','TMEM11','TOMM34','LETM2','LACTB','PNPO','BCKDK',
                 'NT5DC3','SLC25A1','BNIP3','PDK1']
    core8 = ['TMEM11','SLC25A21','GPT2','CYP27A1','MRPS30','ACAA2','GSR','LIG3']
    mitodeg24u = [x.upper() for x in mitodeg24]
    core8u = [x.upper() for x in core8]

    # ---------- 1) threshold-based official vendor results ----------
    thresh_12 = math.log2(1.2)
    thresh_15 = 0.58

    official_dap = vol[vol["Significant"].isin(["Up", "Down"])].copy()
    official_dap["In_MitoCarta3.0"] = official_dap["GeneUpper"].isin(mito_symbols)
    official_dap["In_core8"] = official_dap["GeneUpper"].isin(core8u)
    official_dap["In_MitoDEG24"] = official_dap["GeneUpper"].isin(mitodeg24u)

    strict61 = vol[(vol["FDR"] < 0.05) & (vol["Log2FC"].abs() >= thresh_15)].copy()
    strict61["In_MitoCarta3.0"] = strict61["GeneUpper"].isin(mito_symbols)

    # Gene-centric background and official DAP set
    gene_level_bg = vol.drop_duplicates("GeneUpper")
    bg_total = gene_level_bg["GeneUpper"].nunique()
    bg_mito = int(gene_level_bg["GeneUpper"].isin(mito_symbols).sum())
    dap_gene_set = set(official_dap["GeneUpper"].dropna())
    dap_total = len(dap_gene_set)
    dap_mito = sum(g in mito_symbols for g in dap_gene_set)
    OR, fisher_p = fisher_exact([[dap_mito, dap_total-dap_mito],
                                 [bg_mito-dap_mito, (bg_total-bg_mito)-(dap_total-dap_mito)]],
                                alternative="greater")

    summary = pd.DataFrame({
        "metric": [
            "Vendor-tested proteins in volcano table",
            "Official DAPs at FDR<0.05 and |log2FC|>log2(1.2)",
            "Official DAPs: up",
            "Official DAPs: down",
            "Stricter DAPs at FDR<0.05 and |log2FC|>=0.58",
            "Stricter DAPs: up",
            "Stricter DAPs: down",
            "Gene-level tested background",
            "MitoCarta genes in tested background",
            "Gene-level DAP set",
            "MitoCarta genes in official DAP set",
            "Fisher OR for MitoCarta enrichment",
            "Fisher P for MitoCarta enrichment"
        ],
        "value": [
            len(vol),
            len(official_dap),
            int((official_dap["Significant"]=="Up").sum()),
            int((official_dap["Significant"]=="Down").sum()),
            len(strict61),
            int((strict61["Log2FC"]>0).sum()),
            int((strict61["Log2FC"]<0).sum()),
            bg_total,
            bg_mito,
            dap_total,
            dap_mito,
            OR,
            fisher_p
        ]
    })

    # ---------- 2) global-shift / centering sensitivity ----------
    K_cols = ["K_83","K_85","K_88"]
    W_cols = ["W_80","W_84","W_86"]
    # same testing rule as vendor background: >=2 non-missing per group
    tested_full = prot[(prot[K_cols].notna().sum(axis=1) >= 2) & (prot[W_cols].notna().sum(axis=1) >= 2)].copy()
    sample_log2 = np.log2(tested_full[K_cols + W_cols].astype(float))
    sample_qc = pd.DataFrame({
        "sample": K_cols + W_cols,
        "group": ["KI"]*3 + ["WT"]*3,
        "median_log2_abundance": sample_log2.median(axis=0).values,
        "mean_log2_abundance": sample_log2.mean(axis=0).values
    })

    complete_case = prot[prot[K_cols + W_cols].notna().all(axis=1)].copy()
    complete_case["GeneSymbol"] = complete_case["Gene Name"].apply(first_symbol)
    complete_case["GeneUpper"] = complete_case["GeneSymbol"].astype(str).str.upper()
    logmat = np.log2(complete_case[K_cols + W_cols].astype(float))
    centered = logmat - logmat.median(axis=0)

    complete_case["raw_log2FC"] = logmat[K_cols].mean(axis=1) - logmat[W_cols].mean(axis=1)
    complete_case["centered_log2FC"] = centered[K_cols].mean(axis=1) - centered[W_cols].mean(axis=1)
    pvals = []
    for i in range(centered.shape[0]):
        ki = centered.iloc[i][K_cols].values.astype(float)
        wt = centered.iloc[i][W_cols].values.astype(float)
        pvals.append(ttest_ind(ki, wt, equal_var=True).pvalue)
    complete_case["centered_adj.P.Val"] = multipletests(pvals, method="fdr_bh")[1]
    complete_case["GeneUpper"] = complete_case["GeneUpper"].replace("NAN", np.nan)

    cc_gene = complete_case.dropna(subset=["GeneUpper"]).sort_values("GeneUpper").drop_duplicates("GeneUpper")
    cc_gene["In_MitoCarta3.0"] = cc_gene["GeneUpper"].isin(mito_symbols)
    mito_vals = cc_gene.loc[cc_gene["In_MitoCarta3.0"], "centered_log2FC"]
    non_vals = cc_gene.loc[~cc_gene["In_MitoCarta3.0"], "centered_log2FC"]
    mw_center = mannwhitneyu(mito_vals, non_vals, alternative="greater")

    centered_summary = pd.DataFrame({
        "metric": [
            "Complete-case proteins",
            "Raw median log2FC (KI-WT)",
            "Raw positive fraction",
            "Centered median log2FC (KI-WT)",
            "Centered positive fraction",
            "Proteins with centered adj.P.Val < 0.05",
            "Proteins with centered adj.P.Val < 0.05 and |centered log2FC| > log2(1.2)",
            "Proteins with centered adj.P.Val < 0.05 and |centered log2FC| >= 0.58",
            "Centered median log2FC for MitoCarta proteins",
            "Centered median log2FC for non-MitoCarta proteins",
            "Mann-Whitney P (MitoCarta > non-MitoCarta after centering)"
        ],
        "value": [
            len(complete_case),
            float(np.median(complete_case["raw_log2FC"])),
            float((complete_case["raw_log2FC"] > 0).mean()),
            float(np.median(complete_case["centered_log2FC"])),
            float((complete_case["centered_log2FC"] > 0).mean()),
            int((complete_case["centered_adj.P.Val"] < 0.05).sum()),
            int(((complete_case["centered_adj.P.Val"] < 0.05) & (complete_case["centered_log2FC"].abs() > thresh_12)).sum()),
            int(((complete_case["centered_adj.P.Val"] < 0.05) & (complete_case["centered_log2FC"].abs() >= thresh_15)).sum()),
            float(np.median(mito_vals)),
            float(np.median(non_vals)),
            float(mw_center.pvalue),
        ]
    })

    # ---------- 3) threshold-free rank-based MitoCarta / MitoPathway enrichment ----------
    # signed z statistic from the official vendor differential table
    eps = np.nextafter(0, 1)
    vol["P_value_clipped"] = vol["P value"].clip(lower=eps)
    vol["signed_z"] = np.sign(vol["Log2FC"]) * norm.isf(vol["P_value_clipped"] / 2)
    vol["abs_signed_z"] = vol["signed_z"].abs()
    gene_rank = vol.sort_values(["GeneUpper", "abs_signed_z"], ascending=[True, False]).drop_duplicates("GeneUpper").copy()
    gene_rank["In_MitoCarta3.0"] = gene_rank["GeneUpper"].isin(mito_symbols)

    mito_rank = gene_rank.loc[gene_rank["In_MitoCarta3.0"], "signed_z"]
    non_rank = gene_rank.loc[~gene_rank["In_MitoCarta3.0"], "signed_z"]
    mw_greater = mannwhitneyu(mito_rank, non_rank, alternative="greater")
    auc = mw_greater.statistic / (len(mito_rank) * len(non_rank))

    rank_mitocarta = pd.DataFrame([{
        "n_gene_mapped_tested_proteins": len(gene_rank),
        "n_MitoCarta_genes": int(gene_rank["In_MitoCarta3.0"].sum()),
        "n_non_MitoCarta_genes": int((~gene_rank["In_MitoCarta3.0"]).sum()),
        "rank_metric": "signed_z = sign(log2FC) * norm.isf(P/2)",
        "median_signed_z_MitoCarta": float(np.median(mito_rank)),
        "median_signed_z_non_MitoCarta": float(np.median(non_rank)),
        "Mann_Whitney_two_sided_P": float(mannwhitneyu(mito_rank, non_rank, alternative="two-sided").pvalue),
        "Mann_Whitney_greater_P": float(mw_greater.pvalue),
        "AUC_rank_effect": float(auc),
    }])

    # MitoPathway rank-based enrichment
    score_map = gene_rank.set_index("GeneUpper")["signed_z"]
    universe = set(score_map.index)
    path_rows = []
    for _, row in mitopath.iterrows():
        genes = sorted(set(row["GeneSet"]) & universe)
        if len(genes) < 5:
            continue
        set_scores = score_map.loc[genes]
        bg_scores = score_map.drop(index=genes)
        mw_two = mannwhitneyu(set_scores, bg_scores, alternative="two-sided")
        mw_gt = mannwhitneyu(set_scores, bg_scores, alternative="greater")
        auc_path = mw_gt.statistic / (len(set_scores) * len(bg_scores))
        direction = "positive" if auc_path > 0.5 else "negative"
        path_rows.append({
            "MitoPathway": row["MitoPathway"],
            "Hierarchy": row["MitoPathways Hierarchy"],
            "n_tested_genes": len(genes),
            "median_signed_z_set": float(np.median(set_scores)),
            "median_signed_z_background": float(np.median(bg_scores)),
            "direction": direction,
            "AUC_rank_effect": float(auc_path),
            "P_two_sided": float(mw_two.pvalue),
        })
    rank_paths = pd.DataFrame(path_rows)
    rank_paths["BH_two_sided"] = multipletests(rank_paths["P_two_sided"], method="fdr_bh")[1]
    rank_paths = rank_paths.sort_values(["BH_two_sided", "AUC_rank_effect"], ascending=[True, False])

    # ---------- 4) prespecified family-wise BH ----------
    core7 = vol[vol["GeneUpper"].isin(core8u)].sort_values("GeneUpper").copy()
    core7["BH_within_core7"] = multipletests(core7["P value"], method="fdr_bh")[1]
    core7 = core7.merge(trans[["GeneUpper", "GeneSymbol", "logFC", "P.Value", "adj.P.Val"]], on="GeneUpper", how="left", suffixes=("_prot","_rna"))
    core7["direction_concordant"] = np.sign(core7["Log2FC"]) == np.sign(core7["logFC"])

    mito18 = vol[vol["GeneUpper"].isin(mitodeg24u)].sort_values("GeneUpper").copy()
    mito18["BH_within_MitoDEG18"] = multipletests(mito18["P value"], method="fdr_bh")[1]
    mito18 = mito18.merge(trans[["GeneUpper", "GeneSymbol", "logFC", "P.Value", "adj.P.Val"]], on="GeneUpper", how="left", suffixes=("_prot","_rna"))
    mito18["direction_concordant"] = np.sign(mito18["Log2FC"]) == np.sign(mito18["logFC"])

    # ---------- 5) effect-direction concordance ----------
    all_shared = gene_rank.merge(trans[["GeneUpper", "logFC"]], on="GeneUpper", how="inner")
    all_shared["direction_concordant"] = np.sign(all_shared["Log2FC"]) == np.sign(all_shared["logFC"])
    all_true = int(all_shared["direction_concordant"].sum())
    all_false = int(len(all_shared) - all_true)

    core_true = int(core7["direction_concordant"].sum())
    core_false = int(len(core7) - core_true)
    mito_true = int(mito18["direction_concordant"].sum())
    mito_false = int(len(mito18) - mito_true)

    concord_summary = pd.DataFrame([
        {
            "family": "Core7 formally testable proteins",
            "n": len(core7),
            "concordant": core_true,
            "discordant": core_false,
            "concordance_rate": core_true / len(core7),
            "exact_binomial_P": binomtest(core_true, len(core7), p=0.5, alternative="greater").pvalue,
            "fisher_OR_vs_all_shared_background": fisher_exact([[core_true, core_false], [all_true, all_false]], alternative="greater")[0],
            "fisher_P_vs_all_shared_background": fisher_exact([[core_true, core_false], [all_true, all_false]], alternative="greater")[1],
        },
        {
            "family": "MitoDEG18 formally testable proteins",
            "n": len(mito18),
            "concordant": mito_true,
            "discordant": mito_false,
            "concordance_rate": mito_true / len(mito18),
            "exact_binomial_P": binomtest(mito_true, len(mito18), p=0.5, alternative="greater").pvalue,
            "fisher_OR_vs_all_shared_background": fisher_exact([[mito_true, mito_false], [all_true, all_false]], alternative="greater")[0],
            "fisher_P_vs_all_shared_background": fisher_exact([[mito_true, mito_false], [all_true, all_false]], alternative="greater")[1],
        },
        {
            "family": "All transcriptome-proteome shared genes",
            "n": len(all_shared),
            "concordant": all_true,
            "discordant": all_false,
            "concordance_rate": all_true / len(all_shared),
            "exact_binomial_P": np.nan,
            "fisher_OR_vs_all_shared_background": np.nan,
            "fisher_P_vs_all_shared_background": np.nan,
        },
    ])

    # Write outputs
    outdir = Path(__file__).resolve().parents[1] / "outputs"
    outdir.mkdir(exist_ok=True)
    summary.to_csv(outdir / "proteomics_threshold_summary.csv", index=False)
    official_dap.to_csv(outdir / "official_523_DAPs_with_MitoCarta_and_family_flags.csv", index=False)
    strict61.to_csv(outdir / "strict61_DAPs.csv", index=False)
    sample_qc.to_csv(outdir / "sample_level_log2_abundance_summary.csv", index=False)
    centered_summary.to_csv(outdir / "centered_global_shift_sensitivity_summary.csv", index=False)
    rank_mitocarta.to_csv(outdir / "rank_based_MitoCarta_summary.csv", index=False)
    rank_paths.to_csv(outdir / "rank_based_MitoPathway_enrichment.csv", index=False)
    core7.to_csv(outdir / "core7_family_BH_results.csv", index=False)
    mito18.to_csv(outdir / "MitoDEG18_family_BH_results.csv", index=False)
    concord_summary.to_csv(outdir / "effect_direction_concordance_summary.csv", index=False)
    all_shared.to_csv(outdir / "all_shared_genes_direction_concordance_background.csv", index=False)

    print("Done. Outputs written to:", outdir)

if __name__ == "__main__":
    main()
