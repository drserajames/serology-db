#!/usr/bin/env python3
"""Load the CSVs into a typed DuckDB database + export Parquet.

Produces:
  out/serology.duckdb          persistent DuckDB (antigen/serum/titer_table/titer + titer_flat view)
  out/parquet/<table>/...      Parquet export (titer partitioned by subtype)
"""
import os, duckdb

HERE = os.path.dirname(os.path.abspath(__file__))
ACMACS_DATA = os.environ.get("ACMACS_DATA", os.path.normpath(os.path.join(HERE, os.pardir, "acmacs-data")))
# Output lives outside this code dir (large, WHO-derived); override with SERO_OUT.
OUT = os.environ.get("SERO_OUT", os.path.join(ACMACS_DATA, "serology-db"))
CSV = os.environ.get("SERO_CSV_DIR", os.path.join(OUT, "csv"))
PQ = os.path.join(OUT, "parquet")
DB = os.path.join(OUT, "serology.duckdb")

if os.path.exists(DB):
    os.remove(DB)
con = duckdb.connect(DB)
rd = lambda name: (f"read_csv('{os.path.join(CSV, name)}', header=true, "
                   f"nullstr='', sample_size=-1)")

con.execute(f"""
CREATE TABLE antigen AS SELECT
  ag_id, subtype, name, virus_type, lineage, location, isolation, year,
  passage, reassortant, TRY_CAST(collection_date AS DATE) AS collection_date,
  'hidb' AS source          -- provenance: 'hidb' | 'ace_recovered' (see below)
FROM {rd('antigen.csv')};

CREATE TABLE serum AS SELECT
  sr_id, subtype, serum_id, name, virus_type, lineage, location, isolation,
  year, passage, species, 'hidb' AS source
FROM {rd('serum.csv')};

CREATE TABLE titer_table AS SELECT
  tab_id, subtype, assay, rbc, lab, virus,
  TRY_CAST(table_date AS DATE) AS table_date
FROM {rd('titer_table.csv')};

CREATE TABLE titer AS SELECT
  tab_id, ag_id, sr_id, titer_raw, titer_kind,
  TRY_CAST(titer_value AS INTEGER) AS titer_value,
  TRY_CAST(log_titer AS DOUBLE) AS log_titer,
  -- log_titer is ae Titer::logged() (face value, log2(value/10)). For summary
  -- stats (GMT/SD) the canonical convention is ae Titer::logged_with_thresholded():
  -- a left-censored <N counts as N/2 (log-1), a right-censored >N as N*2 (log+1).
  -- Use THIS column for averaging so censored readings aren't biased.
  CASE titer_kind
    WHEN 'lt' THEN TRY_CAST(log_titer AS DOUBLE) - 1
    WHEN 'gt' THEN TRY_CAST(log_titer AS DOUBLE) + 1
    ELSE TRY_CAST(log_titer AS DOUBLE)
  END AS log_titer_thresholded,
  'hidb' AS source
FROM {rd('titer.csv')};
""")

# Clipped-cell recovery (build_clipped.py, optional): titer cells hidb5 dropped as
# unresolved, recovered from the canonical source .ace charts and tagged provenance
# 'ace_recovered'. Recovered antigens/sera use content-derived ids ({sub}:ra|rs:hash)
# that can't collide with hidb's positional ids; recovered titers keep the original
# tab_id so they attach to existing titer_table rows. Appended here so every
# downstream view/index/export includes them uniformly.
_rec_ag = os.path.join(CSV, "recovered_antigen.csv")
if os.path.exists(_rec_ag):
    con.execute(f"""
    INSERT INTO antigen SELECT
      ag_id, subtype, name, virus_type, lineage, location, isolation, year,
      passage, reassortant, TRY_CAST(collection_date AS DATE), 'ace_recovered'
    FROM {rd('recovered_antigen.csv')};
    INSERT INTO serum SELECT
      sr_id, subtype, serum_id, name, virus_type, lineage, location, isolation,
      year, passage, species, 'ace_recovered'
    FROM {rd('recovered_serum.csv')};
    INSERT INTO titer SELECT
      tab_id, ag_id, sr_id, titer_raw, titer_kind,
      TRY_CAST(titer_value AS INTEGER), TRY_CAST(log_titer AS DOUBLE),
      CASE titer_kind WHEN 'lt' THEN TRY_CAST(log_titer AS DOUBLE) - 1
                      WHEN 'gt' THEN TRY_CAST(log_titer AS DOUBLE) + 1
                      ELSE TRY_CAST(log_titer AS DOUBLE) END,
      'ace_recovered'
    FROM {rd('recovered_titer.csv')};
    """)
    _n = con.execute("SELECT count(*) FROM titer WHERE source='ace_recovered'").fetchone()[0]
    print(f"  recovered   {_n:>9,} clipped titers merged (source='ace_recovered')")

# Sequences: two independent stages joined on seq_id — match.csv (ag_id->seq_id
# + aa, slow, depends on hidb5+seqdb) and clade.csv (seq_id->clade, near-instant,
# depends on clades.json). Kept separate so a clades.json update needn't re-match.
# The `sequence` table is always created so the views below resolve.
match_csv = os.path.join(CSV, "match.csv")
if os.path.exists(match_csv):
    con.execute(f"""
    CREATE TEMP TABLE _match AS SELECT
      ag_id, subtype, seq_id,
      TRY_CAST(aa_length AS INTEGER) AS aa_length,
      TRY_CAST(nuc_length AS INTEGER) AS nuc_length, aa
    FROM {rd('match.csv')};
    """)
    if os.path.exists(os.path.join(CSV, "clade.csv")):
        con.execute(f"CREATE TEMP TABLE _clade AS "
                    f"SELECT seq_id, clade, clade_path FROM {rd('clade.csv')}")
    else:
        con.execute("CREATE TEMP TABLE _clade(seq_id VARCHAR, clade VARCHAR, "
                    "clade_path VARCHAR)")
    con.execute("""
    CREATE TABLE sequence AS SELECT
      m.ag_id, m.subtype, m.seq_id, c.clade, c.clade_path,
      m.aa_length, m.nuc_length, m.aa
    FROM _match m LEFT JOIN _clade c USING (seq_id);
    """)
else:
    con.execute("""CREATE TABLE sequence(
      ag_id VARCHAR, subtype VARCHAR, seq_id VARCHAR, clade VARCHAR,
      clade_path VARCHAR, aa_length INTEGER, nuc_length INTEGER, aa VARCHAR)""")

# Locations resolved via locationdb (build_locations.py). Keyed by the raw hidb
# `O` value in antigen/serum.location. Always created so the views resolve.
if os.path.exists(os.path.join(CSV, "location.csv")):
    con.execute(f"""
    CREATE TABLE location AS SELECT
      location, canonical, country, continent, division,
      TRY_CAST(latitude AS DOUBLE) AS latitude,
      TRY_CAST(longitude AS DOUBLE) AS longitude
    FROM {rd('location.csv')};
    """)
else:
    con.execute("""CREATE TABLE location(
      location VARCHAR, canonical VARCHAR, country VARCHAR, continent VARCHAR,
      division VARCHAR, latitude DOUBLE, longitude DOUBLE)""")

# Convenience flat view (the "one big table" join, materialised on demand).
# Clade is LEFT-joined so titers without a matched sequence still appear.
con.execute("""
CREATE VIEW titer_flat AS
SELECT t.tab_id, tb.subtype, tb.assay, tb.rbc, tb.lab, tb.virus, tb.table_date,
       a.name AS antigen, a.lineage AS antigen_lineage, a.passage AS antigen_passage,
       a.year AS antigen_year, sq.clade AS antigen_clade,
       loc.country AS antigen_country, loc.continent AS antigen_continent,
       s.name AS serum, s.serum_id, s.species AS serum_species,
       t.titer_raw, t.titer_kind, t.titer_value, t.log_titer, t.log_titer_thresholded,
       t.source, a.source AS antigen_source
FROM titer t
JOIN titer_table tb USING (tab_id)
JOIN antigen a USING (ag_id)
JOIN serum  s USING (sr_id)
LEFT JOIN sequence sq USING (ag_id)
LEFT JOIN location loc ON a.location = loc.location;
""")

# Antigen enriched with its sequence/clade + resolved geography (LEFT JOINs keep
# antigens that lack a sequence or an unresolved location).
con.execute("""
CREATE VIEW antigen_sequence AS
SELECT a.*, sq.seq_id, sq.clade, sq.clade_path, sq.aa_length, sq.aa,
       loc.country, loc.continent, loc.latitude, loc.longitude
FROM antigen a
LEFT JOIN sequence sq USING (ag_id)
LEFT JOIN location loc ON a.location = loc.location;
""")

# Helpful indexes for the common lookup paths
for stmt in [
    "CREATE INDEX i_titer_ag ON titer(ag_id)",
    "CREATE INDEX i_titer_sr ON titer(sr_id)",
    "CREATE INDEX i_titer_tab ON titer(tab_id)",
    "CREATE INDEX i_ag_name ON antigen(name)",
    "CREATE INDEX i_ag_loc ON antigen(location)",
    "CREATE INDEX i_sr_name ON serum(name)",
    "CREATE INDEX i_seq_ag ON sequence(ag_id)",
    "CREATE INDEX i_seq_clade ON sequence(clade)",
    "CREATE INDEX i_loc ON location(location)",
    "CREATE INDEX i_loc_country ON location(country)",
]:
    con.execute(stmt)

# Parquet export
os.makedirs(PQ, exist_ok=True)
con.execute(f"COPY antigen     TO '{PQ}/antigen.parquet'     (FORMAT PARQUET)")
con.execute(f"COPY serum       TO '{PQ}/serum.parquet'       (FORMAT PARQUET)")
con.execute(f"COPY titer_table TO '{PQ}/titer_table.parquet' (FORMAT PARQUET)")
con.execute(f"COPY sequence    TO '{PQ}/sequence.parquet'    (FORMAT PARQUET)")
con.execute(f"COPY location    TO '{PQ}/location.parquet'    (FORMAT PARQUET)")
con.execute(f"COPY titer TO '{PQ}/titer' "
            f"(FORMAT PARQUET, PARTITION_BY (subtype), OVERWRITE_OR_IGNORE)"
            if False else
            # titer has no subtype column; join-free partition via tab prefix:
            f"COPY (SELECT *, split_part(tab_id, ':', 1) AS subtype FROM titer) "
            f"TO '{PQ}/titer' (FORMAT PARQUET, PARTITION_BY (subtype), "
            f"OVERWRITE_OR_IGNORE)")

for tbl in ("antigen", "serum", "titer_table", "titer", "sequence", "location"):
    n = con.execute(f"SELECT count(*) FROM {tbl}").fetchone()[0]
    print(f"  {tbl:12s} {n:>9,} rows")
con.close()
print(f"DB: {DB}")
