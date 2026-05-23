"""
PHIND Phage Engineering Dashboard
=================================

Interactive Streamlit dashboard for exploring Phase 1 results: 25 ranked
lytic Listeria phage candidates for reporter-phage engineering.

Run locally:
    streamlit run dashboard/app.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# ---------------------------------------------------------------------------
# Project paths + brand styling
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RANKED_CSV = PROJECT_ROOT / "results" / "candidate_phages_ranked.csv"
CLASSIFIED_CSV = PROJECT_ROOT / "data" / "listeria_phages_classified.csv"
INSERTION_CSV = PROJECT_ROOT / "results" / "insertion_sites_top10.csv"

# PHIND brand palette (matches pitch deck)
TEAL = "#5BA8A0"
TEAL_DARK = "#2D6E68"
NAVY = "#1A2F4D"
LIGHT_BG = "#E8F4F2"

st.set_page_config(
    page_title="PHIND Phage Engineering Atlas",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for brand polish
st.markdown(
    f"""
    <style>
      .block-container {{ padding-top: 2rem; padding-bottom: 2rem; }}
      h1 {{ color: {NAVY}; }}
      h2 {{ color: {TEAL_DARK}; border-bottom: 2px solid {TEAL}; padding-bottom: 4px; }}
      h3 {{ color: {TEAL_DARK}; }}
      .stMetric {{ background-color: {LIGHT_BG}; padding: 12px; border-radius: 8px; border-left: 4px solid {TEAL}; }}
      [data-testid="stMetricLabel"] {{ color: {TEAL_DARK}; font-weight: 600; }}
      [data-testid="stMetricValue"] {{ color: {NAVY}; }}
      div[data-testid="stSidebar"] {{ background-color: {LIGHT_BG}; }}
    </style>
    """,
    unsafe_allow_html=True,
)


# ---------------------------------------------------------------------------
# Data loading (cached so we don't re-read CSV on every interaction)
# ---------------------------------------------------------------------------

@st.cache_data
def load_ranked() -> pd.DataFrame:
    df = pd.read_csv(RANKED_CSV)
    return df


@st.cache_data
def load_full() -> pd.DataFrame:
    df = pd.read_csv(CLASSIFIED_CSV)
    return df


@st.cache_data
def load_insertions() -> pd.DataFrame:
    if INSERTION_CSV.exists():
        return pd.read_csv(INSERTION_CSV)
    return pd.DataFrame()


ranked = load_ranked()
full = load_full()
insertions = load_insertions()


# ---------------------------------------------------------------------------
# Sidebar — project context + filters
# ---------------------------------------------------------------------------

with st.sidebar:
    st.markdown(f"### 🧬 PHIND Phage Atlas")
    st.markdown(
        "An evidence-based ranking of publicly sequenced *Listeria* "
        "phages for **reporter-phage engineering** in the PHIND food-safety "
        "screening platform."
    )
    st.markdown("---")

    st.markdown("#### Filters")

    min_score = st.slider(
        "Minimum engineering readiness score",
        min_value=0,
        max_value=100,
        value=0,
        step=5,
        help="Filter the candidate table to only show phages above this composite score.",
    )

    size_range = st.slider(
        "Genome size range (kb)",
        min_value=int(ranked["length_bp"].min() / 1000),
        max_value=int(ranked["length_bp"].max() / 1000) + 1,
        value=(
            int(ranked["length_bp"].min() / 1000),
            int(ranked["length_bp"].max() / 1000) + 1,
        ),
    )

    # Family selector (parsed from taxonomy)
    def get_family(tax: str) -> str:
        if pd.isna(tax):
            return "Unknown"
        for token in str(tax).split(";"):
            token = token.strip()
            if "viridae" in token.lower():
                return token
        return "Other Caudoviricetes"

    ranked["family"] = ranked["taxonomy"].apply(get_family)
    families = sorted(ranked["family"].unique().tolist())
    selected_families = st.multiselect(
        "Taxonomic family",
        options=families,
        default=families,
    )

    st.markdown("---")
    st.markdown("#### Project")
    st.markdown("[📂 GitHub repo](https://github.com/crystalzys43/PHIND-phage-engineering)")
    st.markdown("Built by **Crystal Zhao** with **Claude** as pair-programmer.")


# Apply filters
filtered = ranked[
    (ranked["engineering_readiness_score"] >= min_score)
    & (ranked["length_bp"] >= size_range[0] * 1000)
    & (ranked["length_bp"] <= size_range[1] * 1000)
    & (ranked["family"].isin(selected_families))
].copy()


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

st.markdown(f"# 🧬 PHIND Phage Engineering Atlas")
st.markdown(
    "**Which phage should PHIND engineer first as a luciferase reporter for "
    "food-safety screening?** This dashboard shows the full Phase 1 pipeline: "
    "from 230 candidate genomes pulled from NCBI down to 25 ranked, "
    "engineering-ready *Listeria* phages."
)

# Funnel metrics
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("NCBI hits", "230", help="Phages matching the keyword search")
c2.metric("Listeria-host", "128", "−102 off-target", delta_color="off")
c3.metric("Lytic candidates", "32", "−96 temperate", delta_color="off")
c4.metric("After dedup", "25", "−7 RefSeq dupes", delta_color="off")
c5.metric("Top-3 mean score", f"{ranked.head(3)['engineering_readiness_score'].mean():.0f}/100")


# ---------------------------------------------------------------------------
# Tabs — main content
# ---------------------------------------------------------------------------

tab1, tab2, tab3, tab4, tab5 = st.tabs(
    [
        "🏆 Ranked Candidates",
        "📊 Score Breakdown",
        "🧬 Genome Landscape",
        "📋 Phage Detail",
        "🎯 Insertion Sites (Phase 2)",
    ]
)


# === Tab 1: Ranked Candidates table ========================================
with tab1:
    st.markdown("## Ranked candidates")
    st.markdown(
        f"Showing **{len(filtered)} of {len(ranked)}** lytic Listeria phages "
        "passing the active filters."
    )

    display = filtered[
        [
            "rank",
            "engineering_readiness_score",
            "accession",
            "name",
            "length_bp",
            "gc_percent",
            "cds_count",
            "family",
            "literature_citation",
        ]
    ].copy()
    display.columns = [
        "Rank", "Score", "Accession", "Phage", "Size (bp)", "GC%",
        "CDS", "Family", "Literature precedent",
    ]

    st.dataframe(
        display,
        hide_index=True,
        use_container_width=True,
        column_config={
            "Rank": st.column_config.NumberColumn(width="small"),
            "Score": st.column_config.ProgressColumn(
                format="%.1f", min_value=0, max_value=100, width="small"
            ),
            "Size (bp)": st.column_config.NumberColumn(format="%d"),
            "Literature precedent": st.column_config.TextColumn(width="large"),
        },
    )

    # Download
    csv_bytes = display.to_csv(index=False).encode()
    st.download_button(
        "⬇ Download filtered CSV",
        csv_bytes,
        file_name="phind_phage_candidates.csv",
        mime="text/csv",
    )

    st.markdown(
        "**How to read this table:** Score 0–100 is a composite across five "
        "biological dimensions (annotation quality, genome size class, lytic "
        "confidence, taxonomic family, and published reporter-phage precedent). "
        "See the **Score Breakdown** tab to inspect any candidate."
    )


# === Tab 2: Score Breakdown per phage =====================================
with tab2:
    st.markdown("## Score breakdown")
    st.markdown(
        "Pick any candidate to see exactly what drove its overall score. "
        "Each component represents a transparent biological criterion."
    )

    phage_options = filtered.apply(
        lambda r: f"#{r['rank']:2d}  {r['accession']:14s}  {r['name'][:50]}",
        axis=1,
    ).tolist()
    selected_label = st.selectbox("Select a phage", phage_options, index=0 if phage_options else None)

    if selected_label:
        idx = phage_options.index(selected_label)
        row = filtered.iloc[idx]

        # Layout: metrics + chart
        col_left, col_right = st.columns([1, 2])

        with col_left:
            st.metric("Total score", f"{row['engineering_readiness_score']:.1f} / 100")
            st.metric("Rank", f"#{int(row['rank'])} of {len(ranked)}")
            st.metric("Genome size", f"{int(row['length_bp']):,} bp")
            st.metric("GC content", f"{row['gc_percent']}%")
            st.metric("Annotated CDS", f"{int(row['cds_count'])}")
            if row["literature_citation"]:
                st.success(f"📚 **Literature precedent:** {row['literature_citation']}")

        with col_right:
            # Bar chart of score components
            components = [
                ("Annotation quality", row["annotation_quality"], 25),
                ("Size class", row["size_class_score"], 25),
                ("Lytic confidence", row["lytic_confidence"], 20),
                ("Taxonomy", row["taxonomy_score"], 20),
                ("Literature bonus", row["literature_bonus"], 10),
            ]

            fig = go.Figure()
            for name, value, max_val in components:
                fig.add_trace(
                    go.Bar(
                        x=[value],
                        y=[name],
                        orientation="h",
                        marker_color=TEAL,
                        name=name,
                        text=f"{value:.1f} / {max_val}",
                        textposition="outside",
                        showlegend=False,
                    )
                )
                # Ghost bar showing max
                fig.add_trace(
                    go.Bar(
                        x=[max_val - value],
                        y=[name],
                        orientation="h",
                        marker_color="rgba(91,168,160,0.15)",
                        showlegend=False,
                        hoverinfo="skip",
                    )
                )

            fig.update_layout(
                barmode="stack",
                height=350,
                margin=dict(l=10, r=80, t=20, b=10),
                xaxis=dict(range=[0, 30], showgrid=False, zeroline=False),
                yaxis=dict(autorange="reversed"),
                plot_bgcolor="white",
                title_text=f"Score components — {row['accession']}",
                title_font_color=NAVY,
            )
            st.plotly_chart(fig, use_container_width=True)


# === Tab 3: Genome Landscape (distribution plots) ==========================
with tab3:
    st.markdown("## Genome landscape across candidates")
    st.markdown(
        "How are the 25 candidates distributed across the variables that matter "
        "for engineering? Use the patterns here to spot natural sub-groups."
    )

    c1, c2 = st.columns(2)

    with c1:
        # Size histogram with score colormap
        fig = px.scatter(
            ranked,
            x="length_bp",
            y="gc_percent",
            color="engineering_readiness_score",
            size="cds_count",
            hover_data=["accession", "name", "rank"],
            color_continuous_scale=[[0, "#CCCCCC"], [1, TEAL]],
            title="Genome size vs GC content (size = annotated CDS)",
            labels={
                "length_bp": "Genome size (bp)",
                "gc_percent": "GC content (%)",
                "engineering_readiness_score": "Readiness score",
            },
        )
        fig.update_layout(plot_bgcolor="white", title_font_color=NAVY)
        st.plotly_chart(fig, use_container_width=True)

    with c2:
        # Family distribution
        family_counts = ranked["family"].value_counts().reset_index()
        family_counts.columns = ["Family", "Count"]
        fig = px.pie(
            family_counts,
            values="Count",
            names="Family",
            title="Taxonomic distribution",
            color_discrete_sequence=[TEAL, TEAL_DARK, NAVY, "#A0D8D2", "#999999"],
        )
        fig.update_layout(title_font_color=NAVY)
        st.plotly_chart(fig, use_container_width=True)

    # Funnel chart for the full pipeline
    st.markdown("### Pipeline funnel")
    funnel = go.Figure(
        go.Funnel(
            y=["NCBI search hits", "Listeria-host", "Lytic candidates", "After dedup"],
            x=[230, 128, 32, 25],
            textinfo="value+percent initial",
            marker={"color": [TEAL, TEAL_DARK, NAVY, "#0F1F33"]},
        )
    )
    funnel.update_layout(
        height=320, margin=dict(l=10, r=10, t=10, b=10), plot_bgcolor="white"
    )
    st.plotly_chart(funnel, use_container_width=True)


# === Tab 4: Phage Detail Card =============================================
with tab4:
    st.markdown("## Phage detail card")
    st.markdown(
        "Full metadata for any candidate, with links to its NCBI record."
    )

    detail_options = ranked.apply(
        lambda r: f"#{r['rank']:2d}  {r['accession']:14s}  {r['name'][:50]}",
        axis=1,
    ).tolist()
    detail_choice = st.selectbox("Choose phage", detail_options, key="detail_select")
    idx = detail_options.index(detail_choice)
    row = ranked.iloc[idx]

    st.markdown(f"### {row['name']}")
    st.markdown(
        f"**Accession:** [`{row['accession']}`]"
        f"(https://www.ncbi.nlm.nih.gov/nuccore/{row['accession']})"
    )

    c1, c2, c3 = st.columns(3)
    c1.metric("Engineering readiness", f"{row['engineering_readiness_score']:.1f}/100")
    c2.metric("Genome size", f"{int(row['length_bp']):,} bp")
    c3.metric("Annotated CDS", f"{int(row['cds_count'])}")

    st.markdown("#### Biology")
    st.markdown(f"- **Organism:** {row.get('organism', '—')}")
    st.markdown(f"- **Host:** {row.get('host', '—') or '—'}")
    st.markdown(f"- **GC content:** {row.get('gc_percent', '—')}%")
    st.markdown(f"- **Lifestyle:** {row.get('lifestyle', '—')}")
    st.markdown(f"- **Isolation source:** {row.get('isolation_source', '—') or '—'}")
    st.markdown(f"- **Country:** {row.get('country', '—') or '—'}")

    st.markdown("#### Taxonomy")
    st.markdown(f"`{row.get('taxonomy', '—')}`")

    if row.get("literature_citation"):
        st.markdown("#### Published reporter-phage precedent")
        st.success(row["literature_citation"])

    st.markdown("#### Score components")
    components_df = pd.DataFrame(
        [
            ("Annotation quality (0–25)", row["annotation_quality"]),
            ("Size class (0–25)", row["size_class_score"]),
            ("Lytic confidence (0–20)", row["lytic_confidence"]),
            ("Taxonomy (0–20)", row["taxonomy_score"]),
            ("Literature bonus (0–10)", row["literature_bonus"]),
            ("**Total**", row["engineering_readiness_score"]),
        ],
        columns=["Component", "Score"],
    )
    st.dataframe(components_df, hide_index=True, use_container_width=True)


# === Tab 5: Insertion Sites (Phase 2) ====================================
with tab5:
    st.markdown("## Reporter cassette insertion sites")
    st.markdown(
        "For each top-ranked phage, this tab shows the highest-scoring "
        "intergenic regions where a luciferase reporter cassette could be "
        "inserted without disrupting the lytic cycle."
    )

    if insertions.empty:
        st.warning(
            "No insertion-site data found. Run "
            "`python src/insertion/01_find_insertion_sites.py` first."
        )
    else:
        st.markdown(
            f"Analyzed **{insertions['accession'].nunique()} phages**, "
            f"found **{len(insertions)}** candidate insertion sites."
        )

        # Phage selector
        phage_acc = st.selectbox(
            "Choose a phage",
            options=insertions["accession"].unique(),
            format_func=lambda a: f"{a}  —  "
            + insertions[insertions['accession'] == a]['name'].iloc[0][:55],
            key="insertion_phage",
        )

        phage_sites = insertions[insertions["accession"] == phage_acc].copy()
        phage_genome_size = full[full["accession"] == phage_acc]["length_bp"].iloc[0]

        # Top 5 table
        st.markdown(f"### Top 5 sites in `{phage_acc}`")
        top5 = phage_sites.head(5)[
            [
                "insertion_score", "site_start", "site_end", "gap_bp",
                "left_gene", "left_category", "right_gene", "right_category",
                "rationale",
            ]
        ]
        top5.columns = [
            "Score", "Start", "End", "Gap (bp)",
            "Left gene", "Left cat.", "Right gene", "Right cat.",
            "Rationale",
        ]
        st.dataframe(
            top5,
            hide_index=True,
            use_container_width=True,
            column_config={
                "Score": st.column_config.ProgressColumn(
                    format="%.1f", min_value=0, max_value=100, width="small"
                ),
                "Gap (bp)": st.column_config.NumberColumn(format="%d"),
                "Rationale": st.column_config.TextColumn(width="large"),
            },
        )

        # Linear genome map
        st.markdown(f"### Linear genome map ({phage_genome_size:,} bp)")
        fig = go.Figure()

        # Genome backbone
        fig.add_shape(
            type="line", x0=0, x1=phage_genome_size, y0=0, y1=0,
            line=dict(color=NAVY, width=3),
        )

        # All sites as faint markers
        for _, s in phage_sites.iterrows():
            color = TEAL if s["insertion_score"] >= 70 else "#999999"
            mid = (s["site_start"] + s["site_end"]) / 2
            fig.add_trace(
                go.Scatter(
                    x=[mid], y=[0.05], mode="markers",
                    marker=dict(size=8, color=color, opacity=0.5),
                    showlegend=False,
                    hovertext=f"Score {s['insertion_score']:.0f}<br>{s['left_gene']} | {s['right_gene']}",
                    hoverinfo="text",
                )
            )

        # Top 5 with labels
        for i, (_, s) in enumerate(phage_sites.head(5).iterrows()):
            mid = (s["site_start"] + s["site_end"]) / 2
            fig.add_trace(
                go.Scatter(
                    x=[mid], y=[0.5 + i * 0.15], mode="markers+text",
                    marker=dict(size=14, color=TEAL_DARK, symbol="triangle-down"),
                    text=[f"#{i+1} ({s['insertion_score']:.0f})"],
                    textposition="top center",
                    textfont=dict(size=10, color=NAVY),
                    showlegend=False,
                    hovertext=f"Score {s['insertion_score']:.0f}<br>"
                    f"Position {int(mid):,}<br>"
                    f"{s['rationale']}",
                    hoverinfo="text",
                )
            )
            # Connector line from label to genome
            fig.add_shape(
                type="line", x0=mid, x1=mid, y0=0, y1=0.5 + i * 0.15,
                line=dict(color=TEAL, width=1, dash="dot"),
            )

        fig.update_layout(
            height=350,
            margin=dict(l=30, r=30, t=30, b=40),
            xaxis=dict(
                title="Genome position (bp)",
                range=[-1000, phage_genome_size + 1000],
                showgrid=False,
            ),
            yaxis=dict(showticklabels=False, showgrid=False, range=[-0.3, 1.5]),
            plot_bgcolor="white",
        )
        st.plotly_chart(fig, use_container_width=True)

        st.info(
            "**Reading the map:** the navy line is the linear phage genome. "
            "Small dots are all candidate intergenic gaps; teal triangles "
            "are the top-5 scoring insertion sites. Hover any marker for "
            "details. Higher-scored sites combine three biological "
            "properties: comfortable gap size, non-essential flanking "
            "genes, and proximity to the lysis gene cluster."
        )

        # Annotation-quality caveat for A511 / P100
        if phage_acc in {"DQ003638.2", "DQ004855.1"}:
            st.warning(
                "⚠ **Annotation caveat:** This phage's NCBI record uses "
                "minimal gene-product names (gp1, gp2, ...) without "
                "functional descriptions. As a result, our keyword-based "
                "categorizer cannot identify essential vs. lysis vs. "
                "permissive flanks confidently, and scores cap around 55. "
                "Phase 2.2 will re-annotate these high-priority candidates "
                "with **Pharokka** to recover functional categories and "
                "produce more accurate insertion-site scoring."
            )


# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------

st.markdown("---")
st.markdown(
    f"<div style='text-align:center; color:{TEAL_DARK}; padding: 10px;'>"
    "PHIND Phage Engineering Atlas · Phase 1 · "
    "<a href='https://github.com/crystalzys43/PHIND-phage-engineering' style='color:"
    f"{TEAL_DARK}'>GitHub</a>"
    "</div>",
    unsafe_allow_html=True,
)
