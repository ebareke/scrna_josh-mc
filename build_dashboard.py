#!/usr/bin/env python3
"""
Build the data layer (assets/js/data.js) and curate/compress figures for the
scRNA Eftud2 / Sf3b4 / Control E9.5 mouse dashboard.

Author: Eric Bareke
"""
import os, json, shutil, subprocess, csv, math
from collections import defaultdict, OrderedDict

SRC = "/Users/ebareke/Documents/Projects/scRNA_Eftud2-Sf3b4-Ctrl"
DST = "/Users/ebareke/Documents/Projects/scrna_josh-mc"
FIG = os.path.join(DST, "assets", "figures")

DATASETS = ["ALL_CELLS", "NEURAL_CREST_CELLS"]
COMPARISONS = ["Eftud2_vs_Control", "Sf3b4_vs_Control", "Eftud2_vs_Sf3b4"]
PADJ = 0.05
LFC = 0.25   # Seurat default-ish effect threshold for "significant" tally

# ----------------------------------------------------------------------------
def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

def fnum(x):
    try:
        if x is None or x == "" or x == "NA":
            return None
        return float(x)
    except (ValueError, TypeError):
        return None

# ----------------------------------------------------------------------------
# Figure curation: copy (and downscale large) PNGs into assets/figures/<ds>/...
# ----------------------------------------------------------------------------
COPIED = {}
def copy_fig(rel_src, rel_dst, max_w=1600):
    """Copy a figure into the repo, downscaling very large PNGs with sips."""
    s = os.path.join(SRC, rel_src)
    if not os.path.exists(s):
        return None
    d = os.path.join(FIG, rel_dst)
    os.makedirs(os.path.dirname(d), exist_ok=True)
    shutil.copy2(s, d)
    # Downscale if the source is large (> ~700 KB) to keep the site light.
    if os.path.getsize(s) > 700_000:
        try:
            subprocess.run(["sips", "--resampleWidth", str(max_w), d],
                           check=True, capture_output=True)
        except Exception as e:
            print("  sips failed for", rel_dst, e)
    COPIED[rel_dst] = True
    return "assets/figures/" + rel_dst

# ----------------------------------------------------------------------------
def build_dataset(ds):
    base = os.path.join(SRC, ds)
    out = {}

    # ---- composition (per-condition %) -------------------------------------
    comp = read_csv(os.path.join(base, "qc_metrics", "conditions_celltype_composition_metrics.csv"))
    cell_types = [r["cell_type"] for r in comp]
    out["composition"] = {
        "cell_types": cell_types,
        "Control": [fnum(r.get("Control")) or 0 for r in comp],
        "Eftud2":  [fnum(r.get("Eftud2"))  or 0 for r in comp],
        "Sf3b4":   [fnum(r.get("Sf3b4"))   or 0 for r in comp],
    }

    # per-replicate composition
    rep = read_csv(os.path.join(base, "qc_metrics", "samples_celltype_composition_metrics.csv"))
    if rep:
        rep_cols = [c for c in rep[0].keys() if c != "cell_type"]
        out["composition_rep"] = {
            "cell_types": [r["cell_type"] for r in rep],
            "samples": rep_cols,
            "values": {c: [fnum(r.get(c)) or 0 for r in rep] for c in rep_cols},
        }

    # ---- composition statistics --------------------------------------------
    chisq = read_csv(os.path.join(base, "composition_stats", "composition_chisq_test.csv"))
    out["chisq"] = chisq
    fisher = read_csv(os.path.join(base, "composition_stats", "composition_fisher_tests.csv"))
    out["fisher"] = [{
        "comparison": r["comparison"], "cell_type": r["cell_type"],
        "OR": fnum(r["OR"]), "log2OR": fnum(r.get("log2OR")),
        "pvalue": fnum(r["pvalue"]), "padj": fnum(r["padj"]),
        "pct_g1": fnum(r.get("pct_g1")), "pct_g2": fnum(r.get("pct_g2")),
    } for r in fisher]
    clr_om = read_csv(os.path.join(base, "composition_stats", "composition_clr_lme_omnibus.csv"))
    out["clr_omnibus"] = [{
        "cell_type": r["cell_type"], "chi2": fnum(r["LRT_chi2"]),
        "pvalue": fnum(r["pvalue"]), "padj": fnum(r["padj"]),
    } for r in clr_om]
    clr_pw = read_csv(os.path.join(base, "composition_stats", "composition_clr_lme_pairwise.csv"))
    out["clr_pairwise"] = [{
        "comparison": r["comparison"], "cell_type": r["cell_type"],
        "delta_clr": fnum(r.get("delta_clr")), "pvalue": fnum(r["pvalue"]),
        "padj": fnum(r["padj"]),
    } for r in clr_pw]

    # ---- velocity summary ---------------------------------------------------
    vel = read_csv(os.path.join(base, "custom_summary_plots", "08_velocity_condition_summary.csv"))
    out["velocity"] = [{
        "condition": r["Condition"], "speed": fnum(r["Speed"]),
        "coherence": fnum(r["Coherence"]), "latent_time": fnum(r["Latent_Time"]),
    } for r in vel]

    # ---- splicing -----------------------------------------------------------
    spl = read_csv(os.path.join(base, "splicing_ratios", "splicing_ratio_by_celltype_condition.csv"))
    out["splicing"] = [{
        "cell_type": r["cell_type"], "condition": r["condition"],
        "mean": fnum(r["mean_ratio"]), "sd": fnum(r["sd_ratio"]), "n": int(float(r["n"])),
    } for r in spl]
    # global per-condition splicing (n-weighted mean across cell types)
    gw = defaultdict(lambda: [0.0, 0])
    for r in out["splicing"]:
        if r["mean"] is not None:
            gw[r["condition"]][0] += r["mean"] * r["n"]
            gw[r["condition"]][1] += r["n"]
    out["splicing_global"] = {c: (v[0]/v[1] if v[1] else None) for c, v in gw.items()}

    # ---- DE per cell type ---------------------------------------------------
    de = {}
    for comp_name in COMPARISONS:
        cdir = os.path.join(base, "de", comp_name)
        if not os.path.isdir(cdir):
            continue
        for ct in sorted(os.listdir(cdir)):
            res = read_csv(os.path.join(cdir, ct, "DE_results.csv"))
            if not res:
                continue
            up = dn = 0
            rows = []
            for r in res:
                lfc = fnum(r.get("avg_log2FC"))
                padj = fnum(r.get("p_val_adj"))
                if lfc is None or padj is None:
                    continue
                sig = (padj < PADJ) and (abs(lfc) >= LFC)
                if sig and lfc > 0: up += 1
                if sig and lfc < 0: dn += 1
                rows.append({"gene": r.get("gene"), "lfc": lfc, "padj": padj,
                             "p": fnum(r.get("p_val")),
                             "pct1": fnum(r.get("pct.1")), "pct2": fnum(r.get("pct.2"))})
            # top genes by significance (smallest padj), keep direction
            rows_sig = [x for x in rows if x["padj"] is not None and x["padj"] < PADJ]
            rows_sig.sort(key=lambda x: (x["padj"], -abs(x["lfc"])))
            top = rows_sig[:25]
            # volcano sample: cap points for payload size, keep all significant + sample of rest
            vol = [x for x in rows if x["padj"] is not None and x["padj"] < 0.1]
            de.setdefault(comp_name, {})[ct] = {
                "up": up, "dn": dn, "n_tested": len(rows),
                "top": top, "volcano": vol[:1500],
            }
    out["de"] = de

    # ---- pseudobulk ---------------------------------------------------------
    pb = {}
    pbdir = os.path.join(base, "pseudobulk")
    if os.path.isdir(pbdir):
        for comp_name in sorted(os.listdir(pbdir)):
            res = read_csv(os.path.join(pbdir, comp_name, "Pseudobulk_DESeq2_results.csv"))
            if not res:
                continue
            up = dn = 0
            rows = []
            for r in res:
                lfc = fnum(r.get("log2FoldChange"))
                padj = fnum(r.get("padj"))
                if lfc is None:
                    continue
                sig = (padj is not None and padj < PADJ) and abs(lfc) >= 1.0
                if sig and lfc > 0: up += 1
                if sig and lfc < 0: dn += 1
                rows.append({"gene": r.get("gene"), "lfc": lfc, "padj": padj,
                             "p": fnum(r.get("pvalue")), "baseMean": fnum(r.get("baseMean"))})
            rows_sig = [x for x in rows if x["padj"] is not None and x["padj"] < PADJ]
            rows_sig.sort(key=lambda x: (x["padj"], -abs(x["lfc"])))
            pb[comp_name] = {"up": up, "dn": dn, "n_tested": len(rows),
                             "top": rows_sig[:30],
                             "volcano": [x for x in rows if x["padj"] is not None][:3000]}
    out["pseudobulk"] = pb

    # ---- enrichment (top GO terms) -----------------------------------------
    def top_terms(path, k=8):
        rows = read_csv(path)
        out_rows = []
        for r in rows:
            out_rows.append({"id": r.get("ID"), "desc": r.get("Description"),
                             "padj": fnum(r.get("p.adjust")), "count": int(float(r.get("Count", 0) or 0)),
                             "ratio": r.get("GeneRatio")})
        out_rows.sort(key=lambda x: (x["padj"] if x["padj"] is not None else 9e9))
        return out_rows[:k]

    enr = {"de": {}, "pseudo_bulk": {}}
    ede = os.path.join(base, "enrichment", "de")
    if os.path.isdir(ede):
        for comp_name in sorted(os.listdir(ede)):
            for ct in sorted(os.listdir(os.path.join(ede, comp_name))):
                ora = os.path.join(ede, comp_name, ct, "ORA")
                if not os.path.isdir(ora):
                    continue
                d = {}
                for onto in ["GO_BP", "GO_CC", "GO_MF"]:
                    for dirn in ["UP", "DOWN"]:
                        fn = f"{ct}_{comp_name}_ORA_{onto}_{dirn}.csv"
                        p = os.path.join(ora, fn)
                        if os.path.exists(p):
                            d[f"{onto}_{dirn}"] = top_terms(p)
                if d:
                    enr["de"].setdefault(comp_name, {})[ct] = d
    epb = os.path.join(base, "enrichment", "pseudo_bulk")
    if os.path.isdir(epb):
        for comp_name in sorted(os.listdir(epb)):
            ora = os.path.join(epb, comp_name, "ORA")
            if not os.path.isdir(ora):
                continue
            d = {}
            for onto in ["GO_BP", "GO_CC", "GO_MF"]:
                for dirn in ["UP", "DOWN"]:
                    fn = f"Pseudobulk_{comp_name}_ORA_{onto}_{dirn}.csv"
                    p = os.path.join(ora, fn)
                    if os.path.exists(p):
                        d[f"{onto}_{dirn}"] = top_terms(p)
            if d:
                enr["pseudo_bulk"][comp_name] = d
    out["enrichment"] = enr

    # ---- figures (curated) --------------------------------------------------
    figs = {}
    # per-condition UMAP (cell type) + condition overlay
    for cond in ["Control", "Eftud2", "Sf3b4"]:
        figs[f"umap_{cond}"] = copy_fig(
            f"{ds}/per_condition/{cond}_umap_cell_type.png",
            f"{ds.lower()}/umap_{cond}.png")
        figs[f"tsne_{cond}"] = copy_fig(
            f"{ds}/per_condition/{cond}_tsne_cell_type.png",
            f"{ds.lower()}/tsne_{cond}.png")
        figs[f"stacked_{cond}"] = copy_fig(
            f"{ds}/per_condition/{cond}_stacked_composition_label.png",
            f"{ds.lower()}/stacked_{cond}.png")
    # side-by-side conditions umap/tsne
    figs["umap_tsne_side"] = copy_fig(
        f"{ds}/qc_metrics/conditions_umap_tsne_majority_side_by_side.png",
        f"{ds.lower()}/conditions_umap_tsne_side.png", max_w=2400)
    # composition stat figures
    for name in ["composition_OR_heatmap", "composition_clr_lme_forest",
                 "conditions_composition_dotplot", "conditions_composition_ranked_bar",
                 "Eftud2_vs_Control_composition_ranked_bar", "Sf3b4_vs_Control_composition_ranked_bar",
                 "Eftud2_vs_Sf3b4_composition_ranked_bar"]:
        figs[name] = copy_fig(f"{ds}/composition_stats/{name}.png",
                              f"{ds.lower()}/{name}.png")
    # custom summary plots
    for name in ["01_composition_dodged_bars", "04_lollipop_both_mutants",
                 "05_bubble_composition_signgene", "06_velocity_metrics_tile_heatmap",
                 "07_delta_velocity_scatter", "08_velocity_summary_bars",
                 "09_sig_genes_per_celltype", "10_shared_genes_deltavelocity_heatmap"]:
        figs[name] = copy_fig(f"{ds}/custom_summary_plots/{name}.png",
                              f"{ds.lower()}/{name}.png")
    # splicing figures
    figs["splicing_by_condition"] = copy_fig(
        f"{ds}/splicing_ratios/splicing_ratio_by_condition.png",
        f"{ds.lower()}/splicing_by_condition.png")
    figs["splicing_by_celltype"] = copy_fig(
        f"{ds}/splicing_ratios/splicing_ratio_by_celltype.png",
        f"{ds.lower()}/splicing_by_celltype.png")
    # pseudobulk volcano
    for comp_name in pb.keys():
        figs[f"pb_volcano_{comp_name}"] = copy_fig(
            f"{ds}/pseudobulk/{comp_name}/Pseudobulk_Volcano.png",
            f"{ds.lower()}/pb_volcano_{comp_name}.png")
    # per-celltype DE volcanoes
    for comp_name in de.keys():
        for ct in de[comp_name].keys():
            key = f"de_volcano_{comp_name}_{ct}"
            figs[key] = copy_fig(
                f"{ds}/de/{comp_name}/{ct}/Volcano.png",
                f"{ds.lower()}/de_volcano_{comp_name}_{ct}.png")
    # enrichment dotplots (GO_BP up/down) per comparison/celltype
    for comp_name in enr["de"].keys():
        for ct in enr["de"][comp_name].keys():
            for dirn in ["UP", "DOWN"]:
                src = f"{ds}/enrichment/de/{comp_name}/{ct}/ORA/{ct}_{comp_name}_ORA_GO_BP_{dirn}_dotplot.png"
                figs[f"enr_dot_{comp_name}_{ct}_{dirn}"] = copy_fig(
                    src, f"{ds.lower()}/enr_dot_{comp_name}_{ct}_{dirn}.png")
    for comp_name in enr["pseudo_bulk"].keys():
        for dirn in ["UP", "DOWN"]:
            src = f"{ds}/enrichment/pseudo_bulk/{comp_name}/ORA/Pseudobulk_{comp_name}_ORA_GO_BP_{dirn}_dotplot.png"
            figs[f"enr_pb_dot_{comp_name}_{dirn}"] = copy_fig(
                src, f"{ds.lower()}/enr_pb_dot_{comp_name}_{dirn}.png")

    out["figures"] = {k: v for k, v in figs.items() if v}
    return out

# ----------------------------------------------------------------------------
def main():
    data = {"datasets": {}, "comparisons": COMPARISONS,
            "conditions": ["Control", "Eftud2", "Sf3b4"],
            "params": {"padj": PADJ, "lfc_de": LFC, "lfc_pb": 1.0}}
    for ds in DATASETS:
        print("Building", ds, "...")
        data["datasets"][ds] = build_dataset(ds)

    js = "// Auto-generated by build_dashboard.py — do not edit by hand.\n"
    js += "window.DASH_DATA = " + json.dumps(data, indent=1) + ";\n"
    out_js = os.path.join(DST, "assets", "js", "data.js")
    with open(out_js, "w") as f:
        f.write(js)
    print("Wrote", out_js, f"({len(js)/1024:.0f} KB)")
    print("Copied", len(COPIED), "figures")

    # also drop the small CSVs into /data for download/provenance
    for ds in DATASETS:
        for rel in ["splicing_ratios/splicing_ratio_by_celltype_condition.csv",
                    "composition_stats/composition_fisher_tests.csv",
                    "composition_stats/composition_clr_lme_pairwise.csv",
                    "composition_stats/composition_clr_lme_omnibus.csv",
                    "qc_metrics/conditions_celltype_composition_metrics.csv",
                    "qc_metrics/samples_celltype_composition_metrics.csv",
                    "custom_summary_plots/08_velocity_condition_summary.csv"]:
            s = os.path.join(SRC, ds, rel)
            if os.path.exists(s):
                d = os.path.join(DST, "data", ds.lower(), rel)
                os.makedirs(os.path.dirname(d), exist_ok=True)
                shutil.copy2(s, d)

if __name__ == "__main__":
    main()
