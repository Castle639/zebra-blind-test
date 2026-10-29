"""Draw the known objects of a blind test from the MPC's own lists, by class, with a seed anyone can change.

    python draw.py MPCORB.DAT.gz CometEls.txt OUT/sample.json [--seed 65] [--since 20261008]
                   [--counts NEO=360,Distant=all,MBA=60,Hilda=20,Trojan=20,other=40] [--comet-years 1.5,2.0] [--t-ref 20261028]
                   [--add 1I,2I,3I]

An asteroid is a candidate when MPCORB's date of last observation (columns 195-202) is on or after --since; its class
comes from MPCORB's a, e (q = a(1-e)): NEO q < 1.3 au; Hilda 3.7 < a < 4.2; Jupiter trojan 5.05 < a < 5.35;
distant a > 5.5 au; main belt 2.0 < a < 3.7; other the rest (Mars-crossers, Hungarias ...). A comet is a candidate when
its perihelion time in CometEls is within the --comet-years window around --t-ref (before, after). Each class is drawn
at random (Python's random.Random(seed), the pool sorted by packed designation first) up to its count; "all" takes the
whole pool. --add names objects taken whatever the pools say (the interstellar objects). Nothing is sent anywhere.
"""
import datetime as dt, gzip, json, random, re, sys

ORDER = ("NEO", "comet", "Distant", "interstellar", "MBA", "Hilda", "Trojan", "other")


def klass(a, e):
    q = a * (1 - e)
    if q < 1.3:
        return "NEO"
    if 3.7 < a < 4.2:
        return "Hilda"
    if 5.05 < a < 5.35:
        return "Trojan"
    if a > 5.5:
        return "Distant"
    if 2.0 < a < 3.7:
        return "MBA"
    return "other"


def comet_tp(l):
    return dt.date(int(l[14:18]), int(l[19:21]), max(1, int(float(l[22:29]))))


def main():
    a = sys.argv
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d  # noqa: E731
    mpcorb, cometels, out = a[1], a[2], a[3]
    seed = int(opt("--seed", "65"))
    since = opt("--since", "20261008")
    counts = dict(kv.split("=") for kv in opt("--counts", "NEO=360,Distant=all,MBA=60,Hilda=20,Trojan=20,other=40").split(","))
    before, after = (float(x) for x in opt("--comet-years", "1.5,2.0").split(","))
    tref = dt.datetime.strptime(opt("--t-ref", "20261028"), "%Y%m%d").date()
    add = [x for x in opt("--add", "1I,2I,3I").split(",") if x]
    rng = random.Random(seed)
    pools = {k: [] for k in counts}
    hdr = True
    with gzip.open(mpcorb, "rt", errors="replace") as fh:
        for l in fh:
            if hdr:
                hdr = not l.startswith("-----")
                continue
            if len(l) < 202 or l[194:202] < since:
                continue
            try:
                A, E, I = float(l[92:103]), float(l[70:79]), float(l[59:68])
            except ValueError:
                continue
            k = klass(A, E)
            if k in pools:
                pools[k].append({"desig": l[166:194].strip(), "packed": l[0:7].strip(), "class": k, "a": A, "e": E,
                                 "q": A * (1 - E), "i": I, "last_obs": l[194:202]})
    sample = []
    for k in counts:
        p = sorted(pools[k], key=lambda x: x["packed"])
        n = len(p) if counts[k] == "all" else min(int(counts[k]), len(p))
        sample += rng.sample(p, n)
    comets = []
    for l in open(cometels):
        if len(l) < 103:
            continue
        try:
            tp = comet_tp(l)
        except ValueError:
            continue
        if re.match(r"^\d+I/", l[102:158].strip()):          # 1I, 2I, 3I: their own class, by --add
            continue
        if -before * 365.25 <= (tp - tref).days <= after * 365.25:
            q, e = float(l[30:39]), float(l[41:49])
            comets.append({"desig": l[102:158].strip(), "packed": l[0:12].strip(), "class": "comet", "q": q, "e": e,
                           "a": q / (1 - e) if e < 1 else None, "i": float(l[71:79]), "tp": tp.isoformat()})
    comets.sort(key=lambda x: x["packed"])
    sample += comets
    for d in add:
        sample.append({"desig": d, "packed": d, "class": "interstellar"})
    for i, o in enumerate(sample):
        o["idx"] = i
    meta = {"seed": seed, "since": since, "counts": counts, "comet_years": [before, after], "t_ref": tref.isoformat(),
            "add": add, "pools": {k: len(v) for k, v in pools.items()}, "comets_in_range": len(comets)}
    json.dump(dict(meta, objects=sample), open(out, "w"), indent=1)
    print("drawn %d: %s; pools %s" % (len(sample), ", ".join("%s %d" % (c, sum(1 for o in sample if o["class"] == c))
                                                            for c in ORDER), meta["pools"]))


if __name__ == "__main__":
    main()
