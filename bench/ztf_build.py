"""Build the blind copy of a ZTF pull (ztf_pull.py): one CSV a linker is given, its truth held out, and its null copy.

    python ztf_build.py PULLDIR OUT/ [--seed 65] [--ss-arcsec 5] [--max-ndethist 2]

The rows: every positive-difference detection that ZTF itself matched to a known object within --ss-arcsec (ztf_ss
ssnamenr, ssdistnr; the truth), and every other positive-difference detection whose ALeRCE object has --max-ndethist or
fewer historical detections (truth unknown). Each row gets a random det_id; no oid, candid, ssnamenr or magnitude-to-name
hint is kept; rows shuffled. Station I41; ALeRCE's ZTF mjd is UTC (the exposure's midpoint, jd - 2400000.5). No tracklet
IDs: zlink pairs the rows itself. The null: the nights' dates shuffled among themselves, each night's rows moved together.

Writes OUT/blind.csv (id, mjd_utc, ra, dec, band, mag, stn), OUT/null.csv, OUT/truth.json (det -> ssnamenr or 0; objects),
OUT/designations.json (each known object's name as Horizons takes it) and OUT/manifest.json.
"""
import csv, glob, hashlib, json, math, os, random, sys

I41_LON = 243.1400


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def night(mjd):
    return int(math.floor(mjd + 2400000.5 + I41_LON / 360.0))


def main():
    a = sys.argv
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d  # noqa: E731
    pull, out = a[1], a[2]
    seed, ssa, maxn = int(opt("--seed", "65")), float(opt("--ss-arcsec", "5")), int(opt("--max-ndethist", "2"))
    os.makedirs(out, exist_ok=True)
    rng = random.Random(seed * 1000 + 41)
    ss = {}
    rd = lambda pat: (r for f in sorted(glob.glob(os.path.join(pull, pat))) for r in csv.DictReader(open(f)))  # noqa: E731
    for r in rd("ss_*.csv"):
        try:
            d = float(r["ssdistnr"])
        except ValueError:
            continue
        nm = (r["ssnamenr"] or "").strip()
        if nm and nm.lower() not in ("null", "none") and 0 <= d <= ssa:
            ss[r["measurement_id"]] = nm
    nd = {r["oid"]: int(float(r["ndethist"] or 0)) for r in rd("obj_*.csv")}
    if not nd or "--oid-count" in a:     # no ztf_object rows pulled: an object's ndethist read as its rows in the pull itself
        from collections import Counter
        nd = Counter(r["oid"] for r in rd("det_*.csv"))
    rows, n_all, n_neg = [], 0, 0
    for r in rd("det_*.csv"):
        n_all += 1
        if str(r["isdiffpos"]).lower() not in ("t", "1", "true"):
            n_neg += 1
            continue
        known = ss.get(r["measurement_id"])
        if known or nd.get(r["oid"], 99) <= maxn:
            mag = r.get("magpsf") or ""
            rows.append((float(r["mjd"]), float(r["ra"]), float(r["dec"]), r["band"], mag, known or 0))
    seen, uniq = set(), []
    for r in rows:
        k = (round(r[0], 7), round(r[1], 7), round(r[2], 7))
        if k not in seen:
            seen.add(k)
            uniq.append(r)
    rows = uniq
    rng.shuffle(rows)
    ids = ["z%07d" % i for i in rng.sample(range(10 ** 7), len(rows))]
    names = sorted({r[5] for r in rows if r[5]})
    code = {nm: i + 1 for i, nm in enumerate(names)}            # the truth as integers, as the Rubin copy's
    truth = {"det": {}, "null_offset": {}, "objects": {}, "names": {str(v): k for k, v in code.items()}}
    hdr = ["id", "mjd_utc", "ra", "dec", "band", "mag", "stn"]
    with open(os.path.join(out, "blind.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(hdr)
        for i, r in zip(ids, rows):
            w.writerow([i, "%.8f" % r[0], "%.7f" % r[1], "%.7f" % r[2], r[3], r[4], "I41"])
            truth["det"][i] = code.get(r[5], 0)
    nights = sorted({night(r[0]) for r in rows})
    perm = list(nights)
    for _ in range(100):
        rng.shuffle(perm)
        if all(p != n for p, n in zip(perm, nights)):
            break
    off = {n: p - n for n, p in zip(nights, perm)}
    truth["null_offset"] = {str(n): off[n] for n in nights}
    with open(os.path.join(out, "null.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(hdr)
        order = list(range(len(rows)))
        rng.shuffle(order)
        for k in order:
            r = rows[k]
            w.writerow([ids[k], "%.8f" % (r[0] + off[night(r[0])]), "%.7f" % r[1], "%.7f" % r[2], r[3], r[4], "I41"])
    per = {}
    for r in rows:
        if r[5]:
            per.setdefault(code[r[5]], {}).setdefault(night(r[0]), 0)
            per[code[r[5]]][night(r[0])] += 1
    for o, d in per.items():
        truth["objects"][str(o)] = {"rows": sum(d.values()), "nights_any": len(d), "nights_2plus": sum(1 for v in d.values() if v >= 2)}
    json.dump(truth, open(os.path.join(out, "truth.json"), "w"))
    # ZTF's ssnamenr: a number ('433') or a provisional designation with no space ('2015XL520'); Horizons takes '2015 XL520'
    def hz(nm):
        if nm.isdigit():
            return nm
        if len(nm) > 6 and nm[:4].isdigit():
            return nm[:4] + " " + nm[4:]
        return nm
    json.dump({str(code[nm]): hz(nm) for nm in names}, open(os.path.join(out, "designations.json"), "w"))
    fo = [o for o, v in truth["objects"].items() if v["nights_2plus"] >= 3]
    man = {"seed": seed, "ss_arcsec": ssa, "max_ndethist": maxn,
           "pull_provenance": [json.loads(l) for l in open(os.path.join(pull, "provenance.jsonl"))],
           "counts": {"det_rows_all": n_all, "negative_difference_left_out": n_neg, "rows": len(rows),
                      "known_rows": sum(1 for r in rows if r[5]), "nights": len(nights), "known_objects": len(per),
                      "findable_2plus_on_3plus_nights": len(fo)},
           "outputs": {f: sha(os.path.join(out, f)) for f in ("blind.csv", "null.csv", "truth.json", "designations.json")}}
    json.dump(man, open(os.path.join(out, "manifest.json"), "w"), indent=1)
    print(json.dumps(man["counts"]))


if __name__ == "__main__":
    main()
