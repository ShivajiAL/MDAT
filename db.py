import sqlite3

DB_NAME = "genotyping.db"

def get_connection():
    return sqlite3.connect(DB_NAME)

def init_db():
    import sqlite3

    con = sqlite3.connect("mabs.db")
    cur = con.cursor()

    # Existing tables (already present in your file)
    cur.execute("""
        CREATE TABLE IF NOT EXISTS uploads (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            label TEXT UNIQUE
        )
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS genotyping (
            upload_id INTEGER,
            marker TEXT,
            line TEXT,
            call TEXT
        )
    """)

    # --------------------------------------------------
    # NEW TABLE: marker position reference (STEP 1)
    # --------------------------------------------------
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
    cur.execute("INSERT OR IGNORE INTO uploads(upload_label) VALUES (?)", (label,))
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
    cur.execute("SELECT upload_id FROM uploads WHERE upload_label=?", (label,))
    res = cur.fetchone()
    con.close()
    return res[0] if res else None
