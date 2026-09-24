import matplotlib
matplotlib.use("Agg")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Patch

priority_list = ["SLE", "eGFR", "UA", "uAlb", "uCr", "uK", "uNa", "SBP"]
trait_color = {
    "SLE": "#d62728",
    "eGFR": "#f67d25",
    "UA": "#9983bd",
    "uAlb": "#5b2d7c",
    "uCr": "#bf6cac",
    "uK": "#89c760",
    "uNa": "#27ae60",
    "SBP": "#fee301"
}

PIE_RADIUS = 0.3

def draw_pie(ax, cx, cy, row_data, all_traits, radius=PIE_RADIUS):
    vals = []
    for tr in all_traits:
        v = row_data.get(tr, 0)
        if pd.isna(v):
            v = 0
        vals.append(int(v))
    total = sum(vals)
    if total <= 0:
        return
    start_angle = 0
    for idx, val in enumerate(vals):
        if val == 0:
            continue
        frac = val / total
        end_angle = start_angle + frac * 360
        wedge = Wedge(
            center=(cx, cy),
            r=radius,
            theta1=start_angle,
            theta2=end_angle,
            facecolor=trait_color[all_traits[idx]],
            edgecolor=None
        )
        ax.add_patch(wedge)
        start_angle = end_angle

#R
meta_df = pd.read_csv("/node1/liuxy/MediCell/trait_compare_2/other_trait/output/plot/meta_info.csv")
plot_wide = pd.read_csv("/node1/liuxy/MediCell/trait_compare_2/other_trait/output/plot/plot_wide.csv")

plot_data = []
for _, row in meta_df.iterrows():
    name = row["name"]
    cell_order = row["cell_order"].split(",")
    gene_list = row["gene_list"].split(",")
    all_traits = row["all_traits"].split(",")
    sub_df = plot_wide[plot_wide["dataset"] == name].copy()
    plot_data.append({
        "cfg": {"name": name, "cell_order": cell_order},
        "plot_wide": sub_df,
        "all_traits": all_traits,
        "gene_list": gene_list
    })

plot_data = [d for d in plot_data if d is not None]
total_rows = len(plot_data)
height_ratios = [1]*total_rows
max_gene_num = max([len(d["gene_list"]) for d in plot_data])
fig_width = max_gene_num * 0.5 * 1.5
fig_height = sum([len(d["cfg"]["cell_order"]) for d in plot_data]) *1.5

fig, axes = plt.subplots(
    nrows=3, ncols=1,
    figsize=(fig_width, fig_height),
    gridspec_kw={"height_ratios": height_ratios, "hspace":0.12},
    squeeze=False
)
axes = axes.flatten()

cell_label_map = {
    "Plasmacytoid_DC": "pDC"
}

for col_idx, data in enumerate(plot_data):
    ax = axes[col_idx]
    cfg = data["cfg"]
    name = cfg["name"]
    cell_order = cfg["cell_order"]
    plot_wide = data["plot_wide"]
    all_traits = data["all_traits"]
    gene_list = data["gene_list"]
    x_max = len(gene_list)
    y_max = len(cell_order)

    for _, row in plot_wide.iterrows():
        draw_pie(ax, row["x"], row["y"], row, all_traits=all_traits)

    ax.set_xlim(0.4, x_max + 0.5)
    ax.set_ylim(y_max + 0.4, 0.6)
    ax.set_xticks(range(1, x_max+1))
    ax.set_xticklabels(gene_list, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(list(range(1, y_max+1)))
    ytick_labels = [cell_label_map.get(ct, ct) for ct in cell_order]
    ax.set_yticklabels(ytick_labels, fontsize=8)
    ax.set_aspect("equal")
    ax.set_title(f"{name} (gene with ≥4 traits in at least one cell)", fontsize=12, pad=6)
    ax.set_xlabel("Shared Genes", fontsize=9)
    if col_idx == 0:
        ax.set_ylabel("Cell Type", fontsize=9)
    else:
        ax.set_ylabel("")

used_traits = set()
for d in plot_data:
    used_traits.update(d["all_traits"])
legend_elements = [Patch(facecolor=trait_color[t], label=t) for t in priority_list if t in used_traits]
fig.legend(
    handles=legend_elements,
    loc="lower center",
    ncol=len(legend_elements),
    bbox_to_anchor=(0.5, 0.04),
    frameon=False,
    fontsize=9
)
plt.tight_layout(rect=[0, 0.08, 1, 0.96])
out_pdf = "/node1/liuxy/MediCell/trait_compare_2/other_trait/output/plot/all_tissue_TWAS_pie_FDR_2.pdf"
plt.savefig(out_pdf, bbox_inches="tight", dpi=300)
plt.close()

