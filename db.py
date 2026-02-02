import sqlite3
import os

# --------------------------------------------------
# DB path: always inside app directory (Cloud-safe)
# --------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(BASE_DIR, "genotyping.db")

def get_connection():
    return sqlite3.connect(DB_PATH)

def init_db():
    con = get_connection()
    cur = con.cursor()

    # Upload metadata
    cur.execute("""
        CREATE TABLE IF NOT EXISTS uploads (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            upload_label TEXT UNIQUE
        )
    """)

    # Genotyping data
    cur.execute("""
        CREATE TABLE IF NOT EXISTS genotyping (
            upload_id INTEGER,
            marker TEXT,
            line TEXT,
            call TEXT
        )
    """)

    # Marker positions (optional, kept as-is)
    cur.execute("""
        CREATE TABLE IF NOT EXISTS marker_positions (
            marker TEXT PRIMARY KEY,
            chr TEXT,
            position_bp INTEGER,
            chr_length_bp INTEGER
        )
    """)

    con.commit()
    con.close()

def add_upload(label):
    con = get_connection()
    cur = con.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO uploads (upload_label) VALUES (?)",
        (label,)
    )
    con.commit()
    con.close()

def get_uploads():
    con = get_connection()
    cur = con.cursor()
    cur.execute("SELECT upload_label FROM uploads")
    res = [r[0] for r in cur.fetchall()]
    con.close()
    return res

def get_upload_id(label):
    con = get_connection()
    cur = con.cursor()
    cur.execute(
        "SELECT id FROM uploads WHERE upload_label=?",
        (label,)
    )
    res = cur.fetchone()
    con.close()
    return res[0] if res else None
