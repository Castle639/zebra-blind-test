"""Build the blind copy of a Rubin pull (rubin_pull.py): one CSV a linker is given, its truth held out, and its null copy.

    python rubin_build.py PULLDIR OUT/ [--seed 65] [--max-sources 4]

The rows: every detection Rubin associated with a known object (sso.csv; truth = its ssObjectId) and every other detection
whose DIA object holds --max-sources sources or fewer (dia.csv with diaobj.csv; truth unknown). Nothing else is kept:
no oid, no measurement_id, no band-to-object hint; each row gets a random det_id and the rows are shuffled. Station X05;
times are ALeRCE's mjd (Rubin's midpoint, TAI), given to zlink with --time-scale tai. No tracklet IDs: the linker pairs
the rows itself (zlink's own pairing, as a tester's file without trkSubs is read).

Writes OUT/blind.csv (det_id, mjd_tai, ra, dec, band), OUT/null.csv (every local night's rows moved together by a whole
number of days, the nights' dates shuffled among themselves: each night keeps its pairs, no object keeps its motion),
OUT/truth.json (det_id -> ssObjectId or 0, the null's offset of each night, the objects' night counts) and OUT/manifest.json.
"""
import csv, hashlib, json, math, os, random, sys

X05_LON = 289.25                       # Rubin's longitude east (ObsCodes: X05), for the local night


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def unpack(p):
    """ssObjectId -> the MPC designation it packs (Rubin packs the MPC's packed designation as 8 ASCII bytes):
    'K00C53G' -> '2000 CG53'; a packed number ('00433', 'A0345', '~0000') -> its number; else None."""
    s = p.to_bytes(8, "big").lstrip(b"\0").decode("ascii", "replace")
    C = {"I": 18, "J": 19, "K": 20}
    if len(s) == 7 and s[0] in C and s[1:3].isdigit():
        cyc = (int(s[4:6]) if s[4].isdigit() else (ord(s[4]) - 55) * 10 + int(s[5])) if s[4:6] != "00" else 0
        return "%d%s %s%s%s" % (C[s[0]], s[1:3], s[3], s[6], cyc or "")
    if len(s) == 5 and s[1:].isdigit():
        return str(int(s) if s[0].isdigit() else (ord(s[0]) - (55 if s[0].isupper() else 61)) * 10000 + int(s[1:]))
    return None


def night(mjd):
    return int(math.floor(mjd + 2400000.5 + X05_LON / 360.0))


def main():
    a = sys.argv
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d  # noqa: E731
    pull, out = a[1], a[2]
    seed, maxs = int(opt("--seed", "65")), int(opt("--max-sources", "4"))
    os.makedirs(out, exist_ok=True)
    rng = random.Random(seed * 1000 + 5)
    small = {r["oid"] for r in csv.DictReader(open(os.path.join(pull, "diaobj.csv"))) if int(float(r["ndiasources"])) <= maxs}
    rows = []
    for r in csv.DictReader(open(os.path.join(pull, "sso.csv"))):
        rows.append((float(r["mjd"]), float(r["ra"]), float(r["dec"]), r["band"], int(r["oid"])))
    n_sso = len(rows)
    n_dia_all = 0
    for r in csv.DictReader(open(os.path.join(pull, "dia.csv"))):
        n_dia_all += 1
        if r["oid"] in small:
            rows.append((float(r["mjd"]), float(r["ra"]), float(r["dec"]), r["band"], 0))
    seen, uniq = set(), []
    for r in rows:
        k = (round(r[0], 7), round(r[1], 7), round(r[2], 7))
        if k not in seen:
            seen.add(k)
            uniq.append(r)
    rows = uniq
    rng.shuffle(rows)
    ids = ["r%07d" % i for i in rng.sample(range(10 ** 7), len(rows))]
    truth = {"det": {}, "null_offset": {}, "objects": {}}
    with open(os.path.join(out, "blind.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "mjd_tai", "ra", "dec", "band", "stn"])
        for i, r in zip(ids, rows):
            w.writerow([i, "%.8f" % r[0], "%.7f" % r[1], "%.7f" % r[2], r[3], "X05"])
            truth["det"][i] = r[4]
    # the null: the nights' dates shuffled among themselves, each night's rows moved together
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
        w.writerow(["id", "mjd_tai", "ra", "dec", "band", "stn"])
        order = list(range(len(rows)))
        rng.shuffle(order)
        for k in order:
            r = rows[k]
            w.writerow([ids[k], "%.8f" % (r[0] + off[night(r[0])]), "%.7f" % r[1], "%.7f" % r[2], r[3], "X05"])
    per = {}
    for r in rows:
        if r[4]:
            per.setdefault(r[4], {}).setdefault(night(r[0]), 0)
            per[r[4]][night(r[0])] += 1
    for o, d in per.items():
        truth["objects"][str(o)] = {"rows": sum(d.values()), "nights_any": len(d), "nights_2plus": sum(1 for v in d.values() if v >= 2)}
    json.dump(truth, open(os.path.join(out, "truth.json"), "w"))
    json.dump({str(o): unpack(o) for o in per if unpack(o)}, open(os.path.join(out, "designations.json"), "w"))
    fo = [o for o, v in truth["objects"].items() if v["nights_2plus"] >= 3]
    man = {"seed": seed, "max_sources": maxs, "pull_provenance": [json.loads(l) for l in open(os.path.join(pull, "provenance.jsonl"))],
           "counts": {"sso_rows": n_sso, "dia_rows_all": n_dia_all, "dia_rows_kept": len(rows) - n_sso, "rows": len(rows),
                      "nights": len(nights), "known_objects": len(per), "findable_2plus_on_3plus_nights": len(fo)},
           "outputs": {f: sha(os.path.join(out, f)) for f in ("blind.csv", "null.csv", "truth.json")}}
    json.dump(man, open(os.path.join(out, "manifest.json"), "w"), indent=1)
    print(json.dumps(man["counts"]))


if __name__ == "__main__":
    main()
