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

        results.append({
            "Plant": plant,
            "RP": F,
            "HET": H,
            "DP": G,
            "NA": I,
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
        "Recovery_%": "BG Recovery %",
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
    auto_hspace = 0.30
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

    plt.tight_layout(rect=[0, 0.06, 1, 0.94])
    return fig
