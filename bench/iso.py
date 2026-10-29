"""The interstellar objects, cut at each night: blind copies with a truth on every row.

    python iso.py ISO.obs80 CROWD.obs80 OUT/ --name 3I --cutoffs 3,4,5,6,8 [--fraction 0.5] [--seed 65] [--window 30]
                  [--obscodes ObsCodes.html]

ISO.obs80: the object's own rows from the MPC (get-obs). CROWD.obs80: the MPC's rows of unnumbered objects of the same
weeks (from MPCAT-OBS's UnnObs.txt.gz), a truth on every row (its designation, columns 6-12). For each cutoff k the window
ends 0.01 d after the last row of the object's k-th local night holding a tracklet of 2+ rows, and holds the object's rows
up to then and a fixed random share (--fraction, by object, seeded) of the crowd's objects. Every row is rewritten as in
build.py (fresh tracklet IDs, the same columns blanked) and shuffled.

Writes OUT/<name>_k<k>.obs80 (among the crowd), OUT/<name>_k<k>_alone.obs80, OUT/<name>_k<k>_null.obs80 (the crowd and the
object, every tracklet moved by its own whole number of days), OUT/<name>_truth.json (every cutoff's t_end, the objects,
the planted keys of every row, the null keys) and OUT/<name>_manifest.json.
"""
import hashlib, json, math, os, random, string, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build  # noqa: E402


def main():
    a = sys.argv
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d  # noqa: E731
    isof, crowdf, out = a[1], a[2], a[3]
    name = opt("--name", "ISO")
    cut = [int(x) for x in opt("--cutoffs", "3,4,5").split(",")]
    frac, seed, window = float(opt("--fraction", "1.0")), int(opt("--seed", "65")), float(opt("--window", "30"))
    lon = build.obscodes(opt("--obscodes", "ObsCodes.html"))
    os.makedirs(out, exist_ok=True)
    rng = random.Random(seed * 1000 + sum(map(ord, name)))
    used = set()

    def fresh():
        while True:
            t = "".join(rng.choice(string.ascii_uppercase + string.digits) for _ in range(7))
            if t not in used and t[0].isalpha():
                used.add(t)
                return t

    def ok(l):
        if len(l) < 80 or l[14] in build.NO_GROUND or l[77:80] not in lon:
            return False
        try:
            build.jd80(l); build.radec(l)
        except ValueError:
            return False
        return True

    night = lambda l: int(math.floor(build.jd80(l) + lon[l[77:80]] / 360.0))  # noqa: E731
    iso = [l.rstrip("\n")[:80] for l in open(isof) if ok(l.rstrip("\n"))]
    g = {}
    for l in iso:
        g.setdefault((l[77:80], night(l)), []).append(l)
    tn = sorted({nt for (s, nt), v in g.items() if len(v) >= 2})
    crowd = [l.rstrip("\n")[:80] for l in open(crowdf) if ok(l.rstrip("\n"))]
    names = sorted({l[0:12] for l in crowd})
    keep = set(rng.sample(names, int(round(frac * len(names)))))
    crowd = [l for l in crowd if l[0:12] in keep]
    isokeys = {build.key(l) for l in iso}
    crowd = [l for l in crowd if build.key(l) not in isokeys]
    objs = [{"idx": 0, "desig": name, "class": "interstellar"}]
    oid = {}
    for nm in sorted({l[0:12] for l in crowd}):
        oid[nm] = len(objs)
        objs.append({"idx": len(objs), "desig": nm.strip(), "class": "crowd"})
    truth = {"name": name, "seed": seed, "fraction": frac, "window": window, "tracklet_nights": tn, "cutoffs": {},
             "objects": objs, "planted_keys": {}, "null_keys": {}}
    # one fresh tracklet ID per object-station-night (the ISO) or designation-station-night (the crowd), shared by all cutoffs
    tid = {}
    rows = []
    for l in iso:
        k = ("ISO", l[77:80], night(l))
        tid.setdefault(k, fresh())
        rows.append((build.norm(l, tid[k]), 0, tid[k]))
    for l in crowd:
        k = (l[0:12], l[77:80], night(l))
        tid.setdefault(k, fresh())
        rows.append((build.norm(l, tid[k]), oid[l[0:12]], tid[k]))
    seen = set()
    rows = [r for r in rows if not (build.key(r[0]) in seen or seen.add(build.key(r[0])))]
    for l, oi, t in rows:
        truth["planted_keys"][build.key(l)] = oi
    man = {"name": name, "iso_sha256": build.sha(isof), "crowd_sha256": build.sha(crowdf), "fraction": frac, "seed": seed,
           "iso_rows": len(iso), "crowd_rows_kept": len(crowd), "crowd_objects_kept": len(oid), "crowd_objects_all": len(names),
           "files": {}}
    for k in cut:
        if k > len(tn):
            continue
        last = max(build.jd80(l) for l in iso if night(l) == tn[k - 1])
        t_end = last + 0.01
        t0 = t_end - window
        win = [r for r in rows if t0 <= build.jd80(r[0]) <= t_end]
        rng.shuffle(win)
        alone = [r for r in win if r[1] == 0]
        # the null: each tracklet its own offset; the ISO's tracklets all different
        by = {}
        for l, oi, t in win:
            by.setdefault(t, []).append((l, oi))
        null = []
        pool = []
        for t in sorted(by):
            jds = [build.jd80(l) for l, _ in by[t]]
            if by[t][0][1] == 0:
                if not pool:
                    pool = list(build.NULL_OFFSETS); rng.shuffle(pool)
                o = next((pool.pop(i) for i in range(len(pool)) if t0 <= min(jds) + pool[i] and max(jds) + pool[i] <= t_end), None)
            else:
                kk, sg = rng.randint(3, 15), rng.choice((1, -1))
                o = next((s * kk for s in (sg, -sg) if t0 <= min(jds) + s * kk and max(jds) + s * kk <= t_end), None)
            if o is None:
                continue
            for l, oi in by[t]:
                nd = max(len(l[23:32].split(".")[1].rstrip()) if "." in l[23:32] else 0, 5)
                nl = l[:15] + build.date80(build.jd80(l) + o, nd)[:17].ljust(17) + l[32:80]
                null.append((nl, oi, t))
                truth["null_keys"][build.key(nl)] = [oi, t, o]
        rng.shuffle(null)
        for suf, rr in (("", win), ("_alone", alone), ("_null", null)):
            fn = "%s_k%d%s.obs80" % (name, k, suf)
            open(os.path.join(out, fn), "w").write("\n".join(r[0] for r in rr) + "\n")
            man["files"][fn] = {"rows": len(rr), "sha256": build.sha(os.path.join(out, fn))}
        ga = {}
        for l, _, t in alone:
            ga.setdefault((t, night(l)), []).append(l)
        iso_nights = sorted({nt for (t, nt), v in ga.items() if len(v) >= 2})
        truth["cutoffs"][str(k)] = {"t_end": t_end, "iso_rows": len(alone), "iso_tracklet_nights": len(iso_nights),
                                    "iso_stations": sorted({l[77:80] for l, _, _ in alone}), "rows": len(win)}
    json.dump(truth, open(os.path.join(out, "%s_truth.json" % name), "w"))
    man["truth_sha256"] = build.sha(os.path.join(out, "%s_truth.json" % name))
    man["cutoffs"] = truth["cutoffs"]
    json.dump(man, open(os.path.join(out, "%s_manifest.json" % name), "w"), indent=1)
    print(json.dumps({k: v for k, v in man.items() if k != "files"}))


if __name__ == "__main__":
    main()
