#!/usr/bin/env python3
"""Resolve antigen/serum locations to country / continent / lat-long via locationdb.

The hidb5 `O` field (our `location` column) is a raw, unresolved mix of full
names ("HONG KONG") and CDC codes ("WI", "AG").  locationdb.json.xz is fully
self-contained, so we resolve directly from it (no ae dependency):

  replacements (spelling)  ->  names (alias)  ->  locations {name: [lat,long,country,division]}
  cdc_abbreviations {code: name}  as a fallback for bare CDC 2-letter codes
  countries {country: continent-index}  ->  continents[]

This matches ae's own name-resolution order (replacements -> names -> locations),
with the CDC-code fallback for the ~13% that name lookup alone misses. ~100%
coverage of antigen+serum locations.

Input : out/csv/{antigen,serum}.csv   (from build_db.py)
Output: out/csv/location.csv          raw location -> country/continent/lat-long
"""
import csv, json, lzma, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ACMACS_DATA = os.environ.get("ACMACS_DATA", "/Users/sarahjames/AC/eu/acmacs-data")
OUT = os.environ.get("SERO_OUT", os.path.join(ACMACS_DATA, "serology-db"))
CSV_DIR = os.environ.get("SERO_CSV_DIR", os.path.join(OUT, "csv"))
LOCDB = os.path.join(ACMACS_DATA, "locationdb.json.xz")


def main():
    if not os.path.exists(LOCDB):
        print(f"!! {LOCDB} missing; skipping locations.", file=sys.stderr)
        return 1
    D = json.loads(lzma.open(LOCDB, "rb").read().decode("utf-8"))
    locations = D["locations"]          # name -> [lat, long, country, division]
    names = D["names"]                  # alias -> canonical name
    replacements = D["replacements"]    # spelling -> name
    cdc = D["cdc_abbreviations"]        # CDC code -> name
    countries = D["countries"]          # country -> continent index
    continents = D["continents"]

    def canon(s):
        if s in replacements:
            s = replacements[s]
        return names.get(s, s)

    def resolve(loc):
        k = canon(loc)
        if k in locations:
            return k, locations[k]
        if loc in cdc:                  # bare CDC 2-letter code fallback
            k2 = canon(cdc[loc])
            if k2 in locations:
                return k2, locations[k2]
        return None, None

    def continent_of(country):
        i = countries.get(country)
        return continents[i] if isinstance(i, int) and 0 <= i < len(continents) else ""

    # distinct locations across antigens + sera
    locs = set()
    for fn in ("antigen.csv", "serum.csv"):
        p = os.path.join(CSV_DIR, fn)
        if not os.path.exists(p):
            continue
        with open(p) as f:
            for r in csv.DictReader(f):
                if r.get("location"):
                    locs.add(r["location"])

    out = open(os.path.join(CSV_DIR, "location.csv"), "w", newline="")
    w = csv.writer(out)
    w.writerow(["location", "canonical", "country", "continent", "division",
                "latitude", "longitude"])
    resolved = 0
    for loc in sorted(locs):
        k, rec = resolve(loc)
        if rec:
            lat, lon, country, division = rec
            w.writerow([loc, k, country, continent_of(country), division, lat, lon])
            resolved += 1
        else:
            w.writerow([loc, "", "", "", "", "", ""])
    out.close()
    print(f"TOTAL  {resolved}/{len(locs)} distinct locations resolved "
          f"({100*resolved/max(len(locs),1):.1f}%)  ->  out/csv/location.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
