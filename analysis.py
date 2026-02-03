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
def analyze_bc(upload_id, rp, dp, bc_samples):
    """
    BC recovery calculation EXACTLY as per A–K logic.
    NA is COUNTED, never derived.
    """

    if not rp or not dp:
        raise ValueError("RP or DP not selected")
    if rp == dp:
        raise ValueError("RP and DP cannot be the same")

    # ----------------------------
    # Load data
    # ----------------------------
    con = get_connection()
    df = pd.read_sql(
        "SELECT marker, line, call FROM genotyping WHERE upload_id=?",
        con,
        params=(upload_id,)
    )
    con.close()

    # ============================
    # PARENT GENOTYPING (A–E)
    # ============================
    parents = df[df.line.isin([rp, dp])].pivot(
        index="marker",
        columns="line",
        values="call"
    )

    A = parents.shape[0]

    parents["status"] = parents.apply(
        lambda r: classify_polymorphism(r[rp], r[dp]),
        axis=1
    )

    na_rp = set(parents[parents[rp] == "NA"].index)
    na_dp = set(parents[parents[dp] == "NA"].index)
    D_set = na_rp.union(na_dp)

    parents_eff = parents.drop(index=D_set)

    poly_markers = set(
    	parents_eff[
        	((parents_eff[rp] == "FAM") & (parents_eff[dp] == "HEX")) |
        	((parents_eff[rp] == "HEX") & (parents_eff[dp] == "FAM"))
    	].index
    )	
    mono_markers = set(parents_eff.index) - poly_markers

    C = len(mono_markers)
    E = A - len(D_set)

    # ============================
    # BC GENOTYPING (F–K)
    # ============================
    bc_df = df[df.line.isin(bc_samples)]
    
    if bc_df.empty:
        raise ValueError(
            "No BC samples found for recovery calculation. "
            "Please upload BC genotyping data first."
        )

    results = []

    for plant, sub in bc_df.groupby("line"):

        # BC calls that ACTUALLY EXIST
        bc_calls = dict(zip(sub.marker, sub.call))

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
            elif bc_call == parents.loc[m, rp]:
                F += 1
            elif bc_call == parents.loc[m, dp]:
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

def get_polymorphic_markers(upload_id, rp, dp):
    con = get_connection()
    df = pd.read_sql(
        "SELECT marker, line, call FROM genotyping WHERE upload_id=?",
        con, params=(upload_id,)
    )
    con.close()

    parents = df[df.line.isin([rp, dp])].pivot(
        index="marker", columns="line", values="call"
    )

    records = []
    for marker, row in parents.iterrows():
        status = classify_polymorphism(row[rp], row[dp])
        if status == "POLYMORPHIC":
            records.append({
                "Marker": marker,
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


def get_marker_status_for_map(upload_id, rp, dp):
    import pandas as pd
    from db import get_connection

    con = get_connection()

    df = pd.read_sql(
        "SELECT marker, line, call FROM genotyping WHERE upload_id=?",
        con,
        params=(upload_id,)
    )

    con.close()

    rp = rp.strip().upper()
    dp = dp.strip().upper()
    df["marker"] = df["marker"].astype(str).str.strip().str.upper()
    df["line"] = df["line"].astype(str).str.strip().str.upper()

    out = []

    for marker, sub in df.groupby("marker"):
        rp_call = sub.loc[sub.line == rp, "call"]
        dp_call = sub.loc[sub.line == dp, "call"]

        status = "OTHER"
        if not rp_call.empty and not dp_call.empty:
            a, b = rp_call.iloc[0], dp_call.iloc[0]
            if (a == "FAM" and b == "HEX") or (a == "HEX" and b == "FAM"):
                status = "POLY"

        out.append({"marker": marker, "status": status})

    return pd.DataFrame(out)


def plot_chromosome_map(marker_status_df, rp, dp):
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
    con = get_connection()
    pos = pd.read_sql("SELECT * FROM marker_positions", con)
    chr_len = pd.read_sql("SELECT * FROM chromosome_lengths", con)
    con.close()

    pos["marker"] = pos["marker"].astype(str).str.upper()
    pos["chr"] = pos["chr"].astype(str)
    chr_len["chr"] = chr_len["chr"].astype(str)
    marker_status_df["marker"] = marker_status_df["marker"].astype(str).str.upper()

    df = marker_status_df.merge(pos, on="marker", how="inner")
    df = df.merge(chr_len, on="chr", how="left")

    df["pos_mb"] = df["position_bp"] / 1e6

    # --------------------------------------------------
    # SAFETY: ensure chromosome length is usable
    # --------------------------------------------------
    if "chr_length_bp" not in df.columns or df["chr_length_bp"].isna().all():
        import streamlit as st
        st.error(
            "Chromosome length information is missing.\n\n"
            "Please upload the chromosome map Excel file "
            "(Sheet 2 must contain: chr, chr_length_bp)."
        )
        return None

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
    auto_hspace = min(base_gap + max_markers_per_chr * density_factor, 1.6)

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
