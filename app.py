import streamlit as st
import pandas as pd
import os
import signal


if "poly_ready" not in st.session_state:
    st.session_state.poly_ready = False

if "show_chr_map" not in st.session_state:
    st.session_state.show_chr_map = False

from db import init_db, get_uploads, get_upload_id
from io_utils import upload_parent_matrix, upload_bc_matrix
from analysis import analyze_bc

# -------------------------------------------------
# Initialize
# -------------------------------------------------
init_db()

st.set_page_config(
    page_title="MABS Tool",
    layout="wide"
)
st.markdown(
    """
    <style>
    /* Limit content width and center it */
    .block-container {
        max-width: 1100px;
        padding-left: 3.5rem;
        padding-right: 3.5rem;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# ---------------------------------------------
# UI CSS (tabs centering + styling)
# ---------------------------------------------
st.markdown(
    """
    <style>
    /* Center align the tab bar */
    div[data-baseweb="tab-list"] {
        justify-content: center;
    }

    /* Increase tab label font size */
    button[data-baseweb="tab"] > div {
        font-size: 18px !important;
        font-weight: 600;
    }

    /* Increase tab padding */
    button[data-baseweb="tab"] {
        padding: 3px 13px !important;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.markdown(
    """
    <style>
    /* Reduce overall vertical spacing */
    .block-container {
        padding-top: 1rem;
        padding-bottom: 0.5rem;
    }

    div[data-testid="stVerticalBlock"] > div {
        gap: 0.4rem;
    }

    h1, h2, h3, h4 {
        margin-bottom: 0.25rem;
    }

    .stSelectbox, .stTextInput, .stFileUploader, .stButton {
        margin-top: 0.25rem;
    }

    /* Disable page scrolling */
    .main {
        overflow-y: hidden;
    }
    </style>
    """,
    unsafe_allow_html=True
)



# -------------------------------------------------
# Header
# -------------------------------------------------
st.markdown(
    "<h1 style='text-align: center;'>🌿 MABS Data Analysis Tool</h1>",
    unsafe_allow_html=True
)


st.markdown("---")

uploads = get_uploads()

# -------------------------------------------------
# MAIN TABS
# -------------------------------------------------
tab1, tab2, tab3, tab4 = st.tabs([
    "🧬 Parent Genotyping data",
    "🔍 Get Polymorphic Markers",
    "🌱 BC Genotyping data",
    "📊 Get BG Recovery & Ranking"
])

# =================================================
# TAB 1: PARENT GENOTYPING
# =================================================
with tab1:

    st.subheader("Upload Parent Genotyping Data")

    label = st.text_input(
        "Upload period (e.g. Jan-2025)",
        key="parent_upload_label"
    )
    file = st.file_uploader(
        "Upload Excel file",
        type="xlsx",
        key="parent_upload_file"
    )

    if st.button("Upload Genotyping Data", key="parent_upload_btn") and file and label:
        upload_parent_matrix(file, label)
        st.success("Parent genotyping data uploaded successfully")
        st.info(f"📅 Upload period: {label}")

    #Block for Marker position data upload

    st.markdown("---")
    st.subheader("Upload SNP Marker Position File")

    pos_file = st.file_uploader(
        "Upload Excel file",
        type="xlsx",
        key="marker_pos_file"
    )


    if st.button("Upload Marker Position Data") and pos_file:
        from analysis import save_marker_position_file

        m, c = save_marker_position_file(pos_file)

        chrom_df = pd.read_excel(pos_file, sheet_name=1)

        # --------------------------------------------------
        # Read chromosome length sheet safely
        # --------------------------------------------------
        chrom_df = pd.read_excel(pos_file, sheet_name=1)

        # Normalize column names
        chrom_df.columns = chrom_df.columns.str.strip().str.lower()

        # Explicit renaming (NO POSITION ASSUMPTIONS)
        rename_map = {}

        for col in chrom_df.columns:
            if col.startswith("chr"):
                rename_map[col] = "chr"
            if "length" in col:
                rename_map[col] = "chr_length_bp"

        chrom_df = chrom_df.rename(columns=rename_map)

        # HARD VALIDATION
        required = {"chr", "chr_length_bp"}
        if not required.issubset(chrom_df.columns):
            st.error(
                f"Chromosome map Sheet 2 must contain columns {required}.\n"
                f"Found columns: {list(chrom_df.columns)}"
            )
            st.stop()

        st.session_state["chrom_df"] = chrom_df
        st.success(f"Marker position data uploaded: {m} markers, {c} chromosomes")


# =================================================
# TAB 2: RP–DP POLYMORPHIC MARKERS
# =================================================
with tab2:

    st.subheader("RP–DP Polymorphic Marker List")

    label = st.selectbox(
        "Genotyping Data Upload Period",
        uploads,
        key="poly_upload_label"
    )

    rp = st.text_input(
        "Recurrent Parent ID",
        key="poly_rp"
    )

    dp = st.text_input(
        "Donor Parent ID",
        key="poly_dp"
    )

    if st.button("View Polymorphic Markers", key="poly_btn"):

        if not rp or not dp:
            st.error("Please enter both RP and DP IDs")
            st.stop()

        uid = get_upload_id(label)

        from analysis import get_polymorphic_markers

        poly = get_polymorphic_markers(uid, rp, dp)

        if poly.empty:
            st.warning("No polymorphic markers found for the selected RP/DP pair")
            st.stop()

        st.metric("🧬 Polymorphic markers identified", len(poly))
        st.dataframe(
            poly,
            height=350,
            use_container_width=True
        )

        # -------------------------------
        # Excel download
        # -------------------------------
        import io

        buffer = io.BytesIO()
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            poly.to_excel(
                writer,
                index=False,
                sheet_name="Polymorphic_Markers"
            )

        buffer.seek(0)

        st.download_button(
            "Download Polymorphic Markers (Excel)",
            buffer,
            file_name="Polymorphic_markers.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

    st.markdown("---")
    
    if st.button("View Chromosome-wise Marker Map"):

        from analysis import get_marker_status_for_map, plot_chromosome_map
   
        
        if "chrom_df" not in st.session_state:
            st.error("Please upload chromosome map file first.")
            st.stop()

        uid = get_upload_id(label)
        
        marker_status = get_marker_status_for_map(uid, rp, dp)

        fig = plot_chromosome_map(
            marker_status,
            rp,
            dp,
            st.session_state["chrom_df"]
        )
        
        st.pyplot(fig)

        import io

        png = io.BytesIO()
        fig.savefig(png, format="png", dpi=300, bbox_inches="tight")
        png.seek(0)

        jpg = io.BytesIO()
        fig.savefig(jpg, format="jpg", dpi=300, bbox_inches="tight")
        jpg.seek(0)

        st.download_button(
            "Download Map (PNG)",
            png,
            file_name=f"Chr_Map_{rp}_{dp}.png",
            mime="image/png"
        )

        st.download_button(
            "Download Map (JPG)",
            jpg,
            file_name=f"Chr_Map_{rp}_{dp}.jpg",
            mime="image/jpeg"
        )


# =================================================
# TAB 3: UPLOAD BC GENOTYPING
# =================================================
with tab3:

    st.subheader("Upload BC Genotyping Data")

    label = st.selectbox(
        "Parent Genotype Upload Period",
        uploads,
        key="bc_upload_label"
    )

    rp = st.text_input(
        "Recurrent Parent ID",
        key="bc_rp"
    )
    dp = st.text_input(
        "Donor Parent ID",
        key="bc_dp"
    )

    file = st.file_uploader(
        "Upload BC Excel file",
        type="xlsx",
        key="bc_upload_file"
    )

    if st.button("Upload Genotyping Data", key="bc_upload_btn") and file:

        if not rp or not dp:
            st.error("Please enter both RP and DP names")
            st.stop()

        bc_df = pd.read_excel(file)
        bc_samples = bc_df.iloc[:, 0].astype(str).tolist()

        st.session_state["bc_samples"] = bc_samples
        

        upload_bc_matrix(file, label)

        st.success(
            f"BC data uploaded ({len(bc_samples)} samples) | RP: {rp} | DP: {dp}"
        )

# =================================================
# TAB 4: BG RECOVERY & RANKING
# =================================================
with tab4:

    # ---------------------------------------------
    # INFO PANEL (OPTION 3 – SECTION PANEL)
    # ---------------------------------------------
    

    st.subheader("BG Recovery & Ranking")

    # ---------------------------------------------
    # INPUTS
    # ---------------------------------------------
    label = st.selectbox(
        "Upload period",
        uploads,
        key="recovery_upload_label"
    )

    rp = st.text_input(
        "Recurrent Parent ID",
        key="recovery_rp"
    )

    dp = st.text_input(
        "Donor Parent ID",
        key="recovery_dp"
    )

    # ---------------------------------------------
    # RUN ANALYSIS
    # ---------------------------------------------
    if st.button("Analyze BG Recovery", key="recovery_btn"):

        uid = get_upload_id(label)
        if "bc_samples" not in st.session_state:
            st.error("Please upload BC genotyping data first")
            st.stop()

        res = analyze_bc(
            uid,
            rp,
            dp,
            st.session_state["bc_samples"]
        )

        # -----------------------------------------
        # OPTION 2: COLORED TABLE
        # -----------------------------------------
        def highlight_top_rank(row):
            if row["Rank"] == 1:
                return ["background-color: #C8E6C9"] * len(row)
            return [""] * len(row)

        styled_res = (
            res.style
            .apply(highlight_top_rank, axis=1)
        )

        st.dataframe(styled_res, use_container_width=True)

        # -----------------------------------------
        # OPTION 4: FORMATTED EXCEL OUTPUT
        # -----------------------------------------
        from io import BytesIO
        from openpyxl import load_workbook
        from openpyxl.styles import Font, PatternFill
        from openpyxl.utils import get_column_letter

        excel_buffer = BytesIO()
        res.to_excel(
            excel_buffer,
            index=False,
            sheet_name="BC_Recovery"
        )
        excel_buffer.seek(0)

        wb = load_workbook(excel_buffer)
        ws = wb.active

        # Bold header
        for cell in ws[1]:
            cell.font = Font(bold=True)

        # Auto column width
        for col in ws.columns:
            max_len = max(
                len(str(cell.value)) if cell.value is not None else 0
                for cell in col
            )
            ws.column_dimensions[
                get_column_letter(col[0].column)
            ].width = max_len + 2

        # Highlight Rank = 1 rows
        green_fill = PatternFill(
            start_color="C8E6C9",
            end_color="C8E6C9",
            fill_type="solid"
        )

        rank_col_idx = None
        for i, cell in enumerate(ws[1], start=1):
            if cell.value == "Rank":
                rank_col_idx = i
                break

        if rank_col_idx:
            for row in ws.iter_rows(min_row=2):
                if row[rank_col_idx - 1].value == 1:
                    for cell in row:
                        cell.fill = green_fill

        final_buffer = BytesIO()
        wb.save(final_buffer)
        final_buffer.seek(0)

        st.download_button(
            label="Download BG Recovery (Excel)",
            data=final_buffer,
            file_name="BC_Recovery_Ranking.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="recovery_download"
        )


# -------------------------------------------------
# Footer
# -------------------------------------------------
st.markdown(
    "<hr style='margin-top:50px;'>"
    "<center><small>Developed by SHIVAJI AJINATH LAVALE</small></center>",
    unsafe_allow_html=True
)


