import pandas as pd
from db import add_upload, get_upload_id, get_connection

VALID_CALLS = ["FAM", "HEX", "HET", "NA"]

def clean_call(x):
    if pd.isna(x):
        return "NA"
    x = str(x).strip().upper()
    return x if x in VALID_CALLS else "NA"

def upload_parent_matrix(file, upload_label):
    df = pd.read_excel(file)
    add_upload(upload_label)
    upload_id = get_upload_id(upload_label)

    marker_col = df.columns[0]
    long = df.melt(id_vars=[marker_col], var_name="line", value_name="call")
    long.rename(columns={marker_col: "marker"}, inplace=True)
    long["call"] = long["call"].apply(clean_call)

    con = get_connection()
    cur = con.cursor()
    for _, r in long.iterrows():
        cur.execute("""
        INSERT OR REPLACE INTO genotyping(upload_id, marker, line, call)
        VALUES (?, ?, ?, ?)
        """, (upload_id, r.marker, r.line, r.call))
    con.commit()
    con.close()

def upload_bc_matrix(file, upload_label):
    df = pd.read_excel(file)
    add_upload(upload_label)
    upload_id = get_upload_id(upload_label)

    plant_col = df.columns[0]
    long = df.melt(id_vars=[plant_col], var_name="marker", value_name="call")
    long.rename(columns={plant_col: "line"}, inplace=True)
    long["call"] = long["call"].apply(clean_call)

    con = get_connection()
    cur = con.cursor()
    for _, r in long.iterrows():
        cur.execute("""
        INSERT OR REPLACE INTO genotyping(upload_id, marker, line, call)
        VALUES (?, ?, ?, ?)
        """, (upload_id, r.marker, r.line, r.call))
    con.commit()
    con.close()
def upload_marker_positions(file):
    import pandas as pd
    import sqlite3

    df = pd.read_excel(file)

    required_cols = {
        "marker",
        "chr",
        "position_bp",
        "chr_length_bp"
    }

    if not required_cols.issubset(df.columns):
        raise ValueError(
            "Excel must contain columns: marker, chr, position_bp, chr_length_bp"
        )

    con = sqlite3.connect("mabs.db")
    cur = con.cursor()

    # Overwrite existing marker position reference
    cur.execute("DELETE FROM marker_positions")

    for _, row in df.iterrows():
        cur.execute(
            """
            INSERT INTO marker_positions
            (marker, chr, position_bp, chr_length_bp)
            VALUES (?, ?, ?, ?)
            """,
            (
                str(row["marker"]),
                str(row["chr"]),
                int(row["position_bp"]),
                int(row["chr_length_bp"])
            )
        )

    con.commit()
    con.close()

    return len(df)

