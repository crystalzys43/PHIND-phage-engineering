"""Generate README figures for the PHIND phage engineering repo."""

from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IMG_DIR = ROOT / "docs" / "images"
IMG_DIR.mkdir(parents=True, exist_ok=True)

TEAL = "#5BA8A0"
TEAL_DARK = "#2D6E68"
NAVY = "#1A2F4D"
LIGHT = "#E8F4F2"
GREY = "#999999"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Arial", "DejaVu Sans"]


# =====================================================================
# Figure 1: Pipeline funnel (labels outside the bars, never clipped)
# =====================================================================
fig, ax = plt.subplots(figsize=(11, 4.2), dpi=150)
fig.patch.set_facecolor("white")

stages = ["NCBI keyword search", "Listeria host (filter 1)",
          "Lytic candidates (filter 2)", "After RefSeq dedup"]
values = [230, 128, 32, 25]
colors = [TEAL, TEAL_DARK, NAVY, "#0F1F33"]

# x-axis grid: bars start at x=0.30 and grow rightward proportionally.
# Labels go on the LEFT (x<0.30), counts on the RIGHT (after the bar),
# loss annotations even further right.
max_val = max(values)
LEFT_BAR = 0.30
BAR_SPAN = 0.45     # bars never exceed 0.45 of figure width
BAR_HEIGHT = 0.55
ROW_GAP = 0.85

for i, (stage, value, color) in enumerate(zip(stages, values, colors)):
    y = -i * ROW_GAP
    width = value / max_val * BAR_SPAN
    rect = mpatches.Rectangle(
        (LEFT_BAR, y), width, BAR_HEIGHT,
        facecolor=color, edgecolor="white", linewidth=2,
    )
    ax.add_patch(rect)

    # Label on the LEFT (right-aligned to bar start)
    ax.text(
        LEFT_BAR - 0.02, y + BAR_HEIGHT / 2, stage,
        ha="right", va="center", fontsize=13, color=NAVY, fontweight="bold",
    )
    # Count immediately AFTER the bar
    ax.text(
        LEFT_BAR + width + 0.015, y + BAR_HEIGHT / 2, f"n = {value}",
        ha="left", va="center", fontsize=13, color=color, fontweight="bold",
    )
    # Loss annotation FURTHER right (only between stages)
    if i < len(stages) - 1:
        loss = values[i] - values[i + 1]
        loss_pct = 100 * loss / values[i]
        ax.annotate(
            f"−{loss}  ({loss_pct:.0f}%)",
            xy=(LEFT_BAR + width / 2, y),
            xytext=(0.95, y - ROW_GAP / 2),
            ha="right", va="center",
            fontsize=10, color="#888888", style="italic",
            arrowprops=dict(arrowstyle="-", color="#cccccc", lw=0.8),
        )

ax.set_xlim(0, 1.0)
ax.set_ylim(-(len(stages) - 1) * ROW_GAP - 0.2, 0.6)
ax.axis("off")
ax.set_title(
    "PHIND Phage Engineering Atlas: 230 candidates  to  25 ranked",
    fontsize=15, color=NAVY, pad=15, fontweight="bold", loc="left", x=0.02,
)
plt.tight_layout()
plt.savefig(IMG_DIR / "pipeline_funnel.png", bbox_inches="tight", facecolor="white")
plt.close()
print("✓ pipeline_funnel.png")


# =====================================================================
# Figure 2: Top-10 scoreboard
# =====================================================================
ranked = pd.read_csv(ROOT / "results" / "candidate_phages_ranked.csv").head(10)

fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
fig.patch.set_facecolor("white")

names = [f"#{int(r['rank'])}  {r['accession']}" for _, r in ranked.iterrows()]
labels = [
    r["name"].replace("Listeria phage ", "").split(",")[0]
    for _, r in ranked.iterrows()
]

components = ["annotation_quality", "size_class_score", "lytic_confidence",
              "taxonomy_score", "literature_bonus"]
comp_labels = ["Annotation\n(0-25)", "Size class\n(0-25)", "Lytic\n(0-20)",
               "Taxonomy\n(0-20)", "Literature\n(0-10)"]
comp_colors = [NAVY, TEAL_DARK, TEAL, "#A0D8D2", "#E6A85C"]

y_positions = list(range(len(ranked) - 1, -1, -1))  # top of chart = rank 1
left = [0] * len(ranked)
for comp, label, color in zip(components, comp_labels, comp_colors):
    ax.barh(
        y_positions, ranked[comp], left=left,
        color=color, label=label, edgecolor="white", height=0.7,
    )
    left = [l + v for l, v in zip(left, ranked[comp])]

# Annotate phage names
for y, label, score in zip(y_positions, labels, ranked["engineering_readiness_score"]):
    ax.text(102, y, f"{label}  ({score:.0f})", va="center",
            fontsize=10, color=NAVY, fontweight="bold")

ax.set_yticks(y_positions)
ax.set_yticklabels(names, fontsize=10, color=NAVY)
ax.set_xlim(0, 130)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("Engineering readiness score (0–100)", fontsize=11, color=NAVY)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(colors=NAVY)

ax.legend(
    loc="lower right",
    bbox_to_anchor=(1.0, -0.18),
    ncol=5, fontsize=9, frameon=False,
)
ax.set_title(
    "Top-10 Listeria phages ranked by reporter-engineering readiness",
    fontsize=14, color=NAVY, pad=15, fontweight="bold",
)
plt.tight_layout()
plt.savefig(IMG_DIR / "top10_scoreboard.png", bbox_inches="tight", facecolor="white")
plt.close()
print("✓ top10_scoreboard.png")


# =====================================================================
# Figure 3: Genome size landscape — colored by lifestyle
# =====================================================================
classified = pd.read_csv(ROOT / "data" / "listeria_phages_classified.csv")

fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
fig.patch.set_facecolor("white")

# Plot temperate
temp = classified[classified["lifestyle"] == "temperate"]
ax.scatter(
    temp["length_bp"] / 1000, temp["gc_percent"],
    s=30, alpha=0.5, color=GREY, label=f"Temperate (n={len(temp)})  — excluded",
    edgecolors="none",
)

# Plot lytic
lytic = classified[classified["lifestyle"] == "lytic"]
ax.scatter(
    lytic["length_bp"] / 1000, lytic["gc_percent"],
    s=60, alpha=0.85, color=TEAL,
    label=f"Lytic (n={len(lytic)})  — reporter candidates",
    edgecolors="white", linewidth=1,
)

# Highlight A511 and P100
for acc, marker_name in [("DQ003638.2", "A511"), ("DQ004855.1", "P100")]:
    row = classified[classified["accession"] == acc]
    if len(row) > 0:
        ax.scatter(
            row["length_bp"].iloc[0] / 1000, row["gc_percent"].iloc[0],
            s=180, marker="*", color="#E6A85C", edgecolors=NAVY, linewidth=1.5,
            zorder=5,
        )
        ax.annotate(
            marker_name,
            (row["length_bp"].iloc[0] / 1000, row["gc_percent"].iloc[0]),
            xytext=(8, 8), textcoords="offset points",
            fontsize=11, fontweight="bold", color=NAVY,
        )

ax.set_xlabel("Genome size (kb)", fontsize=11, color=NAVY)
ax.set_ylabel("GC content (%)", fontsize=11, color=NAVY)
ax.set_title(
    "Listeria phage landscape — lifestyle filter recovers reporter-engineering candidates",
    fontsize=13, color=NAVY, pad=12, fontweight="bold",
)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(colors=NAVY)
ax.grid(True, alpha=0.3)
ax.legend(loc="lower right", fontsize=10, frameon=True, edgecolor="#cccccc")

plt.tight_layout()
plt.savefig(IMG_DIR / "genome_landscape.png", bbox_inches="tight", facecolor="white")
plt.close()
print("✓ genome_landscape.png")

print("\nAll figures saved to", IMG_DIR)
