#!/usr/bin/env python3
"""Demonstrate the query experience over the prototype DB."""
import os, duckdb
_HERE = os.path.dirname(os.path.abspath(__file__))
_ACMACS = os.environ.get("ACMACS_DATA", os.path.normpath(os.path.join(_HERE, os.pardir, "acmacs-data")))
_OUT = os.environ.get("SERO_OUT", os.path.join(_ACMACS, "serology-db"))
con = duckdb.connect(os.path.join(_OUT, "serology.duckdb"), read_only=True)

def show(title, sql):
    print(f"\n### {title}\n" + sql.strip())
    con.sql(sql).show(max_width=110)

show("1. Coverage by subtype x assay (how much data?)", """
SELECT tb.subtype, tb.assay, count(DISTINCT tb.tab_id) AS tables, count(*) AS titers
FROM titer t JOIN titer_table tb USING (tab_id)
GROUP BY ALL ORDER BY titers DESC LIMIT 8
""")

show("2. Most-tested H3 antigens across tables", """
SELECT a.name, a.passage, count(DISTINCT t.tab_id) AS tables, count(*) AS titers
FROM titer t JOIN antigen a USING (ag_id)
WHERE a.subtype='h3'
GROUP BY ALL ORDER BY tables DESC LIMIT 6
""")

show("3. Titer time-series for the most-tested H3 antigen (runtime-chosen)", """
SELECT table_date, lab, assay, titer_raw, log_titer
FROM titer_flat
WHERE subtype='h3' AND titer_kind='num' AND antigen = (
  SELECT antigen FROM titer_flat WHERE subtype='h3' AND titer_kind='num'
  GROUP BY antigen ORDER BY count(*) DESC LIMIT 1)
ORDER BY table_date LIMIT 8
""")

show("4. Geometric-mean reactivity per serum, one table (log->GMT)", """
-- GMT uses log_titer_thresholded (ae's logged_with_thresholded convention):
-- censored readings are included as N/2 (<) or N*2 (>), not dropped.
SELECT serum, serum_species, count(*) n,
       round(pow(2, avg(log_titer_thresholded)) * 10, 1) AS gmt
FROM titer_flat
WHERE tab_id = (SELECT tab_id FROM titer_table WHERE subtype='h3'
               ORDER BY table_date DESC LIMIT 1)
      AND titer_kind IN ('num','lt','gt')
GROUP BY ALL ORDER BY gmt DESC LIMIT 6
""")

show("7. Geographic reactivity: H3 GMT by antigen continent (titers ⨝ locationdb)", """
SELECT antigen_continent AS continent, count(DISTINCT antigen) AS antigens,
       count(*) AS titers, round(pow(2, avg(log_titer_thresholded)) * 10, 0) AS gmt
FROM titer_flat
WHERE subtype='h3' AND antigen_continent IS NOT NULL
      AND titer_kind IN ('num','lt','gt')
GROUP BY ALL ORDER BY titers DESC LIMIT 8
""")

show("6. Clade-aware reactivity: H3 GMT by antigen clade (titers ⨝ seqdb)", """
SELECT antigen_clade AS clade, count(DISTINCT antigen) AS antigens, count(*) AS titers,
       round(pow(2, avg(log_titer_thresholded)) * 10, 0) AS gmt
FROM titer_flat
WHERE subtype='h3' AND antigen_clade IS NOT NULL
      AND titer_kind IN ('num','lt','gt')
GROUP BY ALL HAVING count(*) > 5000 ORDER BY titers DESC LIMIT 8
""")

show("8. Censoring convention matters: H3 HI GMT three ways (data-quality view)", """
-- Same readings, three treatments of left/right-censored titers:
--   drop      = exclude <N / >N entirely (biases GMT upward — drops low titers)
--   facevalue = log_titer          (ae logged(): <N counted as N)
--   canonical = log_titer_thresholded (ae logged_with_thresholded: <N as N/2)
SELECT
  round(pow(2, avg(log_titer)          FILTER(WHERE titer_kind='num')) * 10, 1) AS gmt_drop,
  round(pow(2, avg(log_titer)          FILTER(WHERE titer_kind IN ('num','lt','gt'))) * 10, 1) AS gmt_facevalue,
  round(pow(2, avg(log_titer_thresholded) FILTER(WHERE titer_kind IN ('num','lt','gt'))) * 10, 1) AS gmt_canonical
FROM titer_flat WHERE subtype='h3' AND assay='HI'
""")

show("5. Threshold/censoring breakdown (data-quality view)", """
WITH c AS (
  SELECT tb.subtype, t.titer_kind, count(*) AS n
  FROM titer t JOIN titer_table tb USING (tab_id)
  GROUP BY ALL)
SELECT subtype, titer_kind, n,
       round(100.0*n/sum(n) OVER (PARTITION BY subtype), 1) AS pct
FROM c ORDER BY subtype, n DESC
""")
show("9. Provenance: hidb vs clipped-cell recovery from source charts", """
SELECT split_part(tab_id,':',1) AS subtype, source, count(*) AS titers,
       count(DISTINCT ag_id) AS antigens
FROM titer GROUP BY ALL ORDER BY subtype, source
""")
con.close()
