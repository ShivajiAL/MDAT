import pandas as pd
from db import get_connection

def classify_polymorphism(rp, dp):
    if rp == "NA" or dp == "NA":
        return "NA"
    if rp == "HET" or dp == "HET":
        return "HET"
    if (rp == "FAM" and dp == "HEX") or (rp == "HEX" and dp == "FAM"):
        return "POLYMORPHIC"
    return "MONOMORPHIC"

def analyze_bc(
    bc_df,
    rp,
    dp,
    effective_markers,
    polymorphic_markers
):
    """
    BC recovery calculation EXACTLY as per A–K logic.
    NA is COUNTED, never derived.
    """

    if not rp or not dp:
        raise ValueError("RP or DP not selected")
    if rp == dp:
        raise ValueError("RP and DP cannot be the same")
    
    # ============================
    # PREPARE BC DATA
    # ============================

    bc = bc_df.copy()

    bc.columns = (
        bc.columns.astype(str)
        .str.strip()
    )
 
    line_col = bc.columns[0]

    bc = bc.set_index(line_col)

    bc.index = (
        bc.index.astype(str)
        .str.strip()
    )

    rp = str(rp).strip()
    dp = str(dp).strip()

    bc = bc.applymap(
        lambda x: str(x).strip().upper()
    )
    
    
    if rp not in bc.index:
        raise ValueError(f"{rp} not found in uploaded file.")

    if dp not in bc.index:
        raise ValueError(f"{dp} not found in uploaded file.")

    rp_calls = bc.loc[rp]
    dp_calls = bc.loc[dp]

    poly_markers = []

    mono_markers = []

    for marker in bc.columns:

        rp_call = str(rp_calls[marker]).strip().upper()
        dp_call = str(dp_calls[marker]).strip().upper()

        if rp_call == "NA" or dp_call == "NA":
            continue

        if (
            (rp_call == "FAM" and dp_call == "HEX")
            or
            (rp_call == "HEX" and dp_call == "FAM")
        ):
            poly_markers.append(marker)
        else:
            mono_markers.append(marker)

    C = effective_markers - polymorphic_markers
    E = effective_markers

    # ============================
    # BC GENOTYPING (F–K)
    # ============================
    
    results = []

    for plant in bc.index:
 
        if plant in [rp, dp]:
            continue

        bc_calls = bc.loc[plant].to_dict()

        F = 0  # RP
        G = 0  # DP
        H = 0  # HET
        I = 0  # NA (explicit only)

        # ---- polymorphic markers ----
        for m in poly_markers:
            if m not in bc_calls:
                continue  # missing marker is IGNORED

            bc_call = bc_calls[m]

            if bc_call == "NA":
                I += 1
            elif bc_call == rp_calls[m]:
                F += 1
            elif bc_call == dp_calls[m]:
                G += 1
            elif bc_call == "HET":
                H += 1

        # ---- monomorphic markers ----
        for m in mono_markers:
            if m not in bc_calls:
                continue  # missing marker is IGNORED

            if bc_calls[m] == "NA":
                I += 1

        J = E - I

        if J > 0:
            K = F + (0.5 * H) + C
            rp_recovery = round((K / J) * 100, 2)
        else:
            K = F + (0.5 * H) + C
            rp_recovery = 0
        
        # -----------------------------------------
        # RP% and DP% based on polymorphic markers
        # -----------------------------------------

        poly_total = len(poly_markers)

        poly_valid = poly_total - I

        if poly_valid > 0:
            rp_percent = round(
                ((F + (0.5 * H)) / poly_valid) * 100,
                2
            )

            dp_percent = round(
                ((G + (0.5 * H)) / poly_valid) * 100,
                2
            )
        else:
            rp_percent = 0
            dp_percent = 0

        results.append({
            "Plant": plant,
            "RP": F,
            "HET": H,
            "DP": G,
            "NA": I,
            "RP_%": rp_percent,
            "DP_%": dp_percent,
            "Recovered": K,
            "Total_markers": J,
            "Recovery_%": rp_recovery
        })

    out = pd.DataFrame(results)

    out = out.sort_values(
        by=["RP", "HET", "DP", "NA"],
        ascending=[False, False, True, True]
    ).reset_index(drop=True)

    # sort by ranking criteria
    out = out.sort_values(
    	by=["RP", "HET", "DP", "NA"],
    	ascending=[False, False, True, True]
    ).reset_index(drop=True)

    # dense ranking: same score = same rank, next group = next rank
    out["Rank"] = (
    	out[["RP", "HET", "DP", "NA"]]
    	.apply(tuple, axis=1)
    	.rank(method="dense", ascending=False)
    	.astype(int)
    )

    final_cols = {
        "Plant": "Plant No",
        "RP": "RP",
        "DP": "DP",
        "HET": "HET",
        "NA": "NA",
        "RP_%": "RP %",
        "DP_%": "DP %",
        "Recovery_%": "Overall BG Recovery %",
        "Rank": "Rank"
    }


    out = out[list(final_cols.keys())].rename(columns=final_cols)


    return out

def get_polymorphic_markers(parent_df, rp, dp):

    df = parent_df.copy()

    rp = rp.strip().upper()
    dp = dp.strip().upper()

    # Normalize column names
    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.upper()
    )

    # Normalize marker names
    marker_col = df.columns[0]
    df[marker_col] = (
        df[marker_col]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    if rp not in df.columns:
        raise ValueError(f"RP '{rp}' not found in uploaded file.")

    if dp not in df.columns:
        raise ValueError(f"DP '{dp}' not found in uploaded file.")

    records = []

    for _, row in df.iterrows():

        status = classify_polymorphism(row[rp], row[dp])

        if status == "POLYMORPHIC":

            records.append({

                "Marker": row[marker_col],

                "RP_call": row[rp],

                "DP_call": row[dp]

            })

    return pd.DataFrame(records)

def save_marker_position_file(excel_file):
    import pandas as pd
    from db import get_connection

    # Read sheets
    markers_df = pd.read_excel(excel_file, sheet_name=0)
    chrom_df = pd.read_excel(excel_file, sheet_name=1)

    # Standardize column names
    markers_df.columns = ["marker", "chr", "position_bp"]
    chrom_df.columns = ["chr", "chr_length_bp"]

    # Normalize
    markers_df["marker"] = markers_df["marker"].astype(str).str.strip().str.upper()
    markers_df["chr"] = markers_df["chr"].astype(str).str.strip().str.upper()
    chrom_df["chr"] = chrom_df["chr"].astype(str).str.strip().str.upper()

    con = get_connection()
    cur = con.cursor()

    # Create tables if not exist
    cur.execute("""
        CREATE TABLE IF NOT EXISTS marker_positions (
            marker TEXT PRIMARY KEY,
            chr TEXT,
            position_bp INTEGER
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS chromosome_lengths (
            chr TEXT PRIMARY KEY,
            chr_length_bp INTEGER
        )
    """)

    # Overwrite data
    cur.execute("DELETE FROM marker_positions")
    cur.execute("DELETE FROM chromosome_lengths")

    markers_df.to_sql("marker_positions", con, if_exists="append", index=False)
    chrom_df.to_sql("chromosome_lengths", con, if_exists="append", index=False)

    con.commit()
    con.close()

    return len(markers_df), len(chrom_df)


def get_marker_status_for_map(parent_df, rp, dp):

    df = parent_df.copy()

    rp = rp.strip().upper()
    dp = dp.strip().upper()

    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.upper()
    )

    marker_col = df.columns[0]

    df[marker_col] = (
        df[marker_col]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    if rp not in df.columns:
        raise ValueError(f"RP '{rp}' not found.")

    if dp not in df.columns:
        raise ValueError(f"DP '{dp}' not found.")

    out = []

    for _, row in df.iterrows():

        status = classify_polymorphism(row[rp], row[dp])

        out.append({

            "marker": row[marker_col],

            "status": "POLY" if status == "POLYMORPHIC" else "MONO"

        })

    return pd.DataFrame(out)


def plot_chromosome_map(
    marker_status_df,
    marker_df,
    chrom_df,
    rp,
    dp
):
    
    import pandas as pd
    import matplotlib.pyplot as plt
    from db import get_connection
    from matplotlib.lines import Line2D
    import math
    import re

    # -----------------------------
    # Helper: natural sort
    # -----------------------------
    def natural_key(s):
        return [int(t) if t.isdigit() else t.lower()
                for t in re.findall(r"\d+|\D+", s)]

    # -----------------------------
    # Load marker + chromosome data
    # -----------------------------
    pos = marker_df.copy()

    chr_len = chrom_df.copy()

    pos["marker"] = pos["marker"].astype(str).str.upper()
    pos["chr"] = pos["chr"].astype(str)
    chr_len["chr"] = chr_len["chr"].astype(str)
    marker_status_df["marker"] = marker_status_df["marker"].astype(str).str.upper()

    df = marker_status_df.merge(pos, on="marker", how="inner")
    df = df.merge(chr_len, on="chr", how="left")

    df["pos_mb"] = df["position_bp"] / 1e6
    df["chr_len_mb"] = df["chr_length_bp"] / 1e6

    # -----------------------------
    # Discover chromosomes dynamically
    # -----------------------------
    chromosomes = sorted(df["chr"].unique(), key=natural_key)
    n_chr = len(chromosomes)

    # -----------------------------
    # Decide layout automatically
    # -----------------------------
    if n_chr <= 13:
        n_rows = 1
    elif n_chr <= 26:
        n_rows = 2
    elif n_chr <= 36:
        n_rows = 3
    else:
        n_rows = math.ceil(math.sqrt(n_chr))

    n_cols = math.ceil(n_chr / n_rows)
    
    # -----------------------------
    # Adjust title spacing by layout
    # -----------------------------
    if n_rows == 1:
        title_y = 0.94
        top_rect = 0.88
    else:
        title_y = 0.97
    top_rect = 0.95

    # -----------------------------
    # Global scale & ticks
    # -----------------------------
    global_max = df["chr_len_mb"].max()
    tick_step = 10
    y_ticks = list(range(0, int(global_max) + tick_step, tick_step))

    # -----------------------------
    # Auto row spacing (density-aware)
    # -----------------------------
    max_markers_per_chr = df.groupby("chr")["marker"].count().max()
    base_gap = -0.1
    density_factor = 0.015
    auto_hspace = 0.35

    # -----------------------------
    # Create figure
    # -----------------------------
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(2 * n_cols, 5.5 * n_rows),
        gridspec_kw={"hspace": auto_hspace}
    )

    axes = axes.flatten()

    # -----------------------------
    # Plot each chromosome
    # -----------------------------
    for idx, chrom in enumerate(chromosomes):
        ax = axes[idx]
        sub = df[df["chr"] == chrom].sort_values("pos_mb")

        chr_len_mb = sub["chr_len_mb"].iloc[0]

        # Chromosome backbone
        ax.plot([0, 0], [0, chr_len_mb], color="black", lw=2)

        # Markers
        for i, (_, r) in enumerate(sub.iterrows()):
            color = "red" if r["status"] == "POLY" else "green"
            base_x = 0.08
            final_x = base_x if i % 2 == 0 else -base_x
            ha = "left" if i % 2 == 0 else "right"

            ax.plot(
                [0, final_x],
                [r["pos_mb"], r["pos_mb"]],
                color=color,
                lw=0.4
            )

            ax.text(
                final_x,
                r["pos_mb"],
                r["marker"],
                fontsize=6,
                va="center",
                ha=ha,
                color=color
            )

        # Chromosome label (attached at top)
        ax.text(
            0,
            global_max * -0.015,
            chrom,
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold"
        )


        # Scale handling (first column of each row)
        if idx % n_cols == 0:
            ax.set_yticks(y_ticks)
            ax.set_ylabel("Mb", fontsize=9)
            ax.tick_params(axis="y", labelsize=8)
            ax.spines["left"].set_visible(True)
        else:
            ax.set_yticks([])
            ax.spines["left"].set_visible(False)

        ax.set_xlim(-0.6, 0.6)
        ax.set_ylim(global_max, 0)
        ax.set_xticks([])
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

    # Hide unused axes
    for j in range(len(chromosomes), len(axes)):
        axes[j].axis("off")

    # -----------------------------
    # Title & legend
    # -----------------------------
    fig.suptitle(
        f"Chromosome-Wise Marker Map for {rp} / {dp}",
        fontsize=15,
        fontweight="bold",
        y=0.999
    )

    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               label='Polymorphic markers',
               markerfacecolor='red', markersize=15),
        Line2D([0], [0], marker='o', color='w',
               label='Monomorphic markers',
               markerfacecolor='green', markersize=15)
    ]

    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.02)
    )

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    return fig

def plot_visual_recovery(
    bc_df,
    rp,
    dp,
    marker_df,
    chrom_df,
    ranked_results
):
    """
    Create graphical representation of background recovery.

    - Plants ordered according to BG Recovery ranking
    - Plant IDs shown on left
    - Chromosome widths proportional to actual chromosome length
    - Marker calls converted to continuous segments using midpoint interpolation
    - Chromosome ends inherit color of terminal markers
    """

    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch
    import re

    # -----------------------------------------
    # Helper: natural chromosome sorting
    # -----------------------------------------
    def natural_key(s):
        return [
            int(t) if t.isdigit() else t.lower()
            for t in re.findall(r"\d+|\D+", str(s))
        ]

    # -----------------------------------------
    # Prepare BC genotype data
    # -----------------------------------------
    bc = bc_df.copy()

    bc.columns = (
        bc.columns.astype(str)
        .str.strip()
    )

    line_col = bc.columns[0]

    bc = bc.set_index(line_col)

    bc.index = (
        bc.index.astype(str)
        .str.strip()
    )

    bc = bc.applymap(
        lambda x: str(x).strip().upper()
    )

    rp = str(rp).strip()
    dp = str(dp).strip()

    # -----------------------------------------
    # Prepare marker position information
    # -----------------------------------------
    pos = marker_df.copy()
    chr_len = chrom_df.copy()

    pos["marker"] = (
        pos["marker"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    pos["chr"] = pos["chr"].astype(str).str.strip()

    chr_len["chr"] = (
        chr_len["chr"]
        .astype(str)
        .str.strip()
    )

    pos["position_bp"] = pd.to_numeric(
        pos["position_bp"],
        errors="coerce"
    )

    chr_len["chr_length_bp"] = pd.to_numeric(
        chr_len["chr_length_bp"],
        errors="coerce"
    )

    # -----------------------------------------
    # Use all markers from SNP position file
    # -----------------------------------------
    pos["marker"] = (
        pos["marker"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    # Merge chromosome lengths
    pos = pos.merge(
        chr_len[["chr", "chr_length_bp"]],
        on="chr",
        how="left"
    )

    # Remove incomplete positional records
    pos = pos.dropna(
        subset=["position_bp", "chr_length_bp"]
    )

    # -----------------------------------------
    # Chromosome order
    # -----------------------------------------
    chromosomes = sorted(
        pos["chr"].unique(),
        key=natural_key
    )

    # -----------------------------------------
    # Get plants in BG recovery rank order
    # -----------------------------------------
    ranked = ranked_results.copy()

    # Your analyze_bc output uses "Plant No"
    # -----------------------------------------
    # Plant order
    # RP first, DP second, followed by BC plants
    # in BG recovery rank order
    # -----------------------------------------
    bc_plant_order = (
        ranked
        .sort_values("Rank")
        ["Plant No"]
        .astype(str)
        .str.strip()
        .tolist()
    )

    plant_order = [rp, dp] + bc_plant_order

    # Keep only plants actually present
    plant_order = [
        p for p in plant_order
        if p in bc.index
    ]

    # -----------------------------------------
    # Classification colors
    # -----------------------------------------
    colors = {
        "RP": "#2E8B57",   # green
        "DP": "#E31A1C",   # red
        "HET": "#FFD92F",  # yellow
        "NA": "#FFFFFF",   # white
        "MONO": "#90EE90"      # light green
    }

    # -----------------------------------------
    # Calculate chromosome width ratios
    # -----------------------------------------
    chromosome_lengths = []

    for chrom in chromosomes:

        length = (
            pos.loc[
                pos["chr"] == chrom,
                "chr_length_bp"
            ]
            .iloc[0]
        )

        chromosome_lengths.append(length)

    # -----------------------------------------
    # Figure sizing
    # -----------------------------------------
    n_plants = len(plant_order)
    
    # -----------------------------------------
    # Row spacing
    # -----------------------------------------
    row_height = 0.85      # Actual coloured row height
    row_spacing = 1.0      # Distance between rows

    y_positions = []

    for i, plant in enumerate(plant_order):

        # Normal spacing
        y = i * row_spacing

        # Add extra separation after RP and DP
        if i >= 2:
            y += 0.40

        y_positions.append(y)
    
    fig_height = max(
        5,
        0.38 * n_plants + 2
    )

    fig_width = 24

    fig = plt.figure(
        figsize=(fig_width, fig_height)
    )

    gs = GridSpec(
        1,
        len(chromosomes),
        figure=fig,
        width_ratios=chromosome_lengths,
        wspace=0.08
    )

    axes = []

    # -----------------------------------------
    # Create chromosome axes
    # -----------------------------------------
    for i, chrom in enumerate(chromosomes):

        ax = fig.add_subplot(gs[0, i])

        axes.append(ax)

        chr_markers = (
            pos[pos["chr"] == chrom]
            .sort_values("position_bp")
            .reset_index(drop=True)
        )

        chr_length = (
            chr_markers["chr_length_bp"]
            .iloc[0]
        )

        marker_names = chr_markers["marker"].tolist()
        marker_positions = chr_markers["position_bp"].tolist()

        # -------------------------------------
        # Plot each plant
        # -------------------------------------
        for row_idx, plant in enumerate(plant_order):

            # Y position for this plant
            y = y_positions[row_idx]
            
            calls = []

            for marker in marker_names:

            # -----------------------------------------
            # Monomorphic marker
            # Marker is not present in BC genotype file
            # -----------------------------------------
                if marker not in bc.columns:

                    status = "MONO"

                else:

                    plant_call = str(
                        bc.loc[plant, marker]
                    ).strip().upper()

                    rp_call = str(
                        bc.loc[rp, marker]
                    ).strip().upper()

                    dp_call = str(
                        bc.loc[dp, marker]
                    ).strip().upper()

                    if plant_call == "NA":

                        status = "NA"

                    elif plant_call == "HET":

                        status = "HET"

                    elif plant_call == rp_call:

                        status = "RP"

                    elif plant_call == dp_call:

                        status = "DP"

                    else:

                        status = "NA"

                calls.append(status)

            # ---------------------------------
            # Midpoint interpolation boundaries
            # ---------------------------------
            boundaries = [0]

            for j in range(len(marker_positions) - 1):
                midpoint = (
                    marker_positions[j]
                    + marker_positions[j + 1]
                ) / 2

                boundaries.append(midpoint)

            boundaries.append(chr_length)

            # ---------------------------------
            # Draw genotype segments
            # ---------------------------------
            for j, status in enumerate(calls):

                start = boundaries[j]
                end = boundaries[j + 1]

                rect = patches.Rectangle(
                    (start, y),
                    end - start,
                    row_height,
                    facecolor=colors[status],
                    edgecolor="none"
                )

                ax.add_patch(rect)

        # -------------------------------------
        # Axis formatting
        # -------------------------------------
        ax.set_xlim(0, chr_length)
        ax.set_ylim(
            y_positions[-1] + row_height,
            -0.1
        )

        ax.set_xticks([])
        ax.set_yticks([])

        ax.set_title(
            chrom,
            fontsize=9,
            fontweight="bold",
            pad=8
        )

        # Remove axis borders
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Plant IDs only on first chromosome
        if i == 0:

            ax.set_yticks(
                [
                    y + row_height / 2
                    for y in y_positions
                ]
            )

            ax.set_yticklabels(
                plant_order,
                fontsize=8
            )

            ax.tick_params(
                axis="y",
                length=0,
                pad=5
            )

    # -----------------------------------------
    # Main title
    # -----------------------------------------
    fig.suptitle(
        f"Visual Representation of Background Recovery ({rp} × {dp})",
        fontsize=15,
        fontweight="bold",
        y=0.98
    )

    # -----------------------------------------
    # Legend
    # -----------------------------------------
    legend_elements = [
        Patch(facecolor=colors["RP"], edgecolor="black", label="RP"),
        Patch(facecolor=colors["DP"], edgecolor="black", label="DP"),
        Patch(facecolor=colors["HET"], edgecolor="black", label="HET"),
        Patch(facecolor=colors["NA"], edgecolor="black", label="NA"),
        Patch(facecolor=colors["MONO"], edgecolor="black", label="Monomorphic")
    ]

    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, 0.01)
    )

    plt.tight_layout(
        rect=[0.08, 0.06, 1, 0.94]
    )

    return fig
