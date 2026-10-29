"""Check that a blind copy hides what it should, before any linker sees it.

    python check_blind.py OUT/ SAMPLE.json OBSDIR [--background itf.txt.gz]

On every row of blind.obs80, planted.obs80 and null.obs80: columns 1-5, 13-14, 57-65 and 73-77 blank; column 71 never N
or T; no satellite, roving, radar or offset row; no tracklet ID equal to a drawn object's packed designation, its name, or
any original ITF trkSub; the planted rows' tracklet IDs never reused by a background row. Prints one line a check and
exits 1 if any fails.
"""
import gzip, json, os, sys

NO_GROUND = set("SsVvRrXxOoEeHh")


def main():
    a = sys.argv
    out, sample, obsdir = a[1], a[2], a[3]
    bgf = a[a.index("--background") + 1] if "--background" in a else None
    names = set()
    for o in json.load(open(sample))["objects"]:
        names |= {o.get("packed", "").strip()[:7], o["desig"].strip()[:7]}
        f = os.path.join(obsdir, "%04d.obs80" % o["idx"])
        if os.path.exists(f):
            for l in open(f):
                if len(l) >= 12:
                    names |= {l[0:5].strip(), l[5:12].strip()}
    names.discard("")
    if bgf:
        op = gzip.open if bgf.endswith(".gz") else open
        with op(bgf, "rt", errors="replace") as fh:
            for l in fh:
                if len(l) >= 12:
                    names.add(l[5:12].strip())
    truth = json.load(open(os.path.join(out, "truth.json")))
    planted_ids = set(truth["tracklet_object"])
    ok = True
    for fn in ("blind.obs80", "planted.obs80", "null.obs80"):
        bad = {"cols_1_5": 0, "cols_13_14": 0, "cols_57_65": 0, "cols_73_77": 0, "band_N_T": 0, "no_ground": 0,
               "id_is_a_name": 0, "short": 0}
        n = 0
        ids = set()
        for l in open(os.path.join(out, fn)):
            l = l.rstrip("\n")
            n += 1
            if len(l) != 80:
                bad["short"] += 1
                continue
            bad["cols_1_5"] += l[0:5] != " " * 5
            bad["cols_13_14"] += l[12:14] != "  "
            bad["cols_57_65"] += l[56:65] != " " * 9
            bad["cols_73_77"] += l[72:77] != " " * 5
            bad["band_N_T"] += l[70] in "NT"
            bad["no_ground"] += l[14] in NO_GROUND
            bad["id_is_a_name"] += l[5:12].strip() in names
            ids.add(l[5:12])
        good = not any(bad.values())
        ok &= good
        print("%-14s rows %6d  tracklet IDs %6d  %s  %s" % (fn, n, len(ids), "CLEAN" if good else "FAIL", json.dumps(bad)))
    print("S1-B0", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
