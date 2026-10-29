"""Build a blind copy: known objects' MPC rows, disguised, mixed into a frozen background, plus its null copy.

    python build.py SAMPLE.json OBSDIR OUT/ --t-end JD [--window 30] [--seed 65] [--background itf.txt.gz] [--obscodes ObsCodes.html]

Every row, planted or background, is written the same way (MPC 80-column):
  columns 1-5 blank; 6-12 a fresh random tracklet ID (one per object-station-night for a planted object, one per
  original trkSub-station-night for a background row); 13-14 blank (discovery mark, note 1); 15 kept (the observation
  type); 16-56 kept (time, RA, Dec); 57-65 blank; 66-70 kept (magnitude); 71 kept unless N or T (a comet's nuclear or
  total magnitude, blanked: it would say "comet"); 72 kept (the star catalogue); 73-77 blank (the publication
  reference); 78-80 kept (the station). Satellite, roving, radar and offset rows are left out (no ground position).
  A planted row equal to a background row in time, place and station is one detection: the planted copy is kept.
  The rows are shuffled.

Writes OUT/:
  blind.obs80      the planted rows among the background rows of the window (what a linker is given)
  blind.csv        the same rows as CSV (det_id = "blind.obs80:<line>", obsTime ISO UTC, jd_utc, ra_deg, dec_deg, mag,
                   band, stn, trkSub), for any linker
  planted.obs80    the planted rows alone (the same tracklet IDs), shuffled
  null.obs80       the null copy: blind.obs80 with every tracklet moved by its own whole number of days (3 to 15, either
                   sign, kept inside the window); a planted object's tracklets get different offsets where it has 26 or
                   fewer of them, so no real object is left in it
  truth.json       held out: the objects (class, elements, rows, tracklets, nights, findable) and for each file the
                   object of every planted row by its time-place-station key, and the null copy's offsets
  manifest.json    the parameters, every source's and output's sha256, counts
Nothing is sent anywhere.
"""
import csv, datetime as dt, gzip, hashlib, json, math, os, random, string, sys

NO_GROUND = set("SsVvRrXxOoEeHh")
NULL_OFFSETS = [s * k for k in range(3, 16) for s in (1, -1)]       # 26 offsets in days


def jd80(l):
    y, m, d = int(l[15:19]), int(l[20:22]), float(l[23:32])
    a = (14 - m) // 12; yy = y + 4800 - a; mm = m + 12 * a - 3
    return int(d) + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045 - 0.5 + (d - int(d))


def date80(jd, ndec):
    """jd -> 'YYYY MM DD.dddddd' with ndec decimals of the day (the row's own precision kept)."""
    n = int(round((jd + 0.5) * 10 ** ndec)); Z, frac = n // 10 ** ndec, n % 10 ** ndec
    a = int((Z - 1867216.25) / 36524.25); A = Z + 1 + a - a // 4
    B = A + 1524; C = int((B - 122.1) / 365.25); D = int(365.25 * C); E = int((B - D) / 30.6001)
    day = B - D - int(30.6001 * E)
    mon = E - 1 if E < 14 else E - 13
    yr = C - 4716 if mon > 2 else C - 4715
    return "%04d %02d %02d.%0*d" % (yr, mon, day, ndec, frac)


def radec(l):
    h, mi, s = float(l[32:34]), float(l[35:37]), float(l[38:44])
    sg = -1.0 if l[44] == "-" else 1.0
    dd, dm, ds = float(l[45:47]), float(l[48:50]), float(l[51:56])
    return 15.0 * (h + mi / 60 + s / 3600), sg * (dd + dm / 60 + ds / 3600)


def obscodes(path):
    lon = {}
    for l in open(path, errors="replace"):
        if len(l) > 30 and l[3] == " ":
            try:
                lon[l[0:3]] = float(l[4:13])
                float(l[13:21]); float(l[21:30])
            except ValueError:
                continue
    return lon


def key(l):
    return l[15:56] + l[77:80]


def norm(l, trk):
    l = l[:80].ljust(80)
    band = l[70] if l[70] not in "NT" else " "
    return " " * 5 + trk + "  " + l[14] + l[15:56] + " " * 9 + l[65:70] + band + l[71] + " " * 5 + l[77:80]


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    a = sys.argv
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d  # noqa: E731
    sample, obsdir, out = a[1], a[2], a[3]
    t_end, window, seed = float(opt("--t-end", None)), float(opt("--window", "30")), int(opt("--seed", "65"))
    bgf = opt("--background", None)
    lon = obscodes(opt("--obscodes", "ObsCodes.html"))
    os.makedirs(out, exist_ok=True)
    rng = random.Random(seed * 1000 + 65)
    t0 = t_end - window
    used = set()
    # a fresh ID once met an old ITF trkSub by chance (R27J7AQ, station 644 in 2007); every ID ever in
    # the background file or in the planted objects' rows is now refused, so no fresh ID is anyone's name
    forbid = set()
    for f in [bgf] if bgf else []:
        op = gzip.open if f.endswith(".gz") else open
        with op(f, "rt", errors="replace") as fh:
            for l in fh:
                forbid.add(l[5:12].strip())
    for o in json.load(open(sample))["objects"]:
        forbid |= {o.get("packed", "").strip()[:7], o["desig"].strip()[:7]}
        f = os.path.join(obsdir, "%04d.obs80" % o["idx"])
        if os.path.exists(f):
            for l in open(f):
                forbid |= {l[0:5].strip(), l[5:12].strip()}

    def fresh():
        while True:
            t = "".join(rng.choice(string.ascii_uppercase + string.digits) for _ in range(7))
            if t not in used and t not in forbid and t[0].isalpha():
                used.add(t)
                return t

    def night(jd, stn):
        return int(math.floor(jd + lon[stn] / 360.0))

    def usable(l):
        if len(l) < 80 or l[14] in NO_GROUND or l[77:80] not in lon:
            return False
        try:
            jd = jd80(l); radec(l)
        except ValueError:
            return False
        return t0 <= jd <= t_end

    # ---- planted rows, one fresh ID an object-station-night
    S = json.load(open(sample))
    objs, planted, pkey = [], [], {}
    for o in S["objects"]:
        f = os.path.join(obsdir, "%04d.obs80" % o["idx"])
        rec = dict(o, fetched=os.path.exists(f), rows=0, tracklets=0, nights=0, tracklets_any=0, findable=False)
        if rec["fetched"]:
            raw = [l.rstrip("\n") for l in open(f) if usable(l.rstrip("\n"))]
            groups = {}
            for l in raw:
                if key(l) in pkey:                                  # a row already planted (an MPC duplicate)
                    continue
                pkey[key(l)] = o["idx"]
                groups.setdefault((l[77:80], night(jd80(l), l[77:80])), []).append(l)
            nts = set()
            for (stn, nt), ls in sorted(groups.items()):
                t = fresh()
                for l in ls:
                    planted.append((norm(l, t), o["idx"], t))
                if len(ls) >= 2:
                    nts.add(nt)
                    rec["tracklets"] += 1
            rec["rows"] = sum(len(v) for v in groups.values())
            rec["tracklets_any"] = len(groups)
            rec["nights"] = len(nts)
            rec["nights_any"] = len({nt for (_, nt) in groups})
            rec["findable"] = rec["tracklets"] >= 3 and len(nts) >= 3
        objs.append(rec)
    pkeys = {key(l) for l, _, _ in planted}

    # ---- background rows of the window, one fresh ID an original trkSub-station-night
    bg, bg_dup, bg_in = [], 0, 0
    if bgf:
        op = gzip.open if bgf.endswith(".gz") else open
        groups = {}
        d_lo, d_hi = date80(t0 - 1, 0)[:10], date80(t_end + 1, 0)[:10]   # 'YYYY MM DD', a fast first cut
        with op(bgf, "rt", errors="replace") as fh:
            for l in fh:
                l = l.rstrip("\n")
                if len(l) < 80 or not (d_lo <= l[15:25] <= d_hi):
                    continue
                if not usable(l):
                    continue
                bg_in += 1
                if key(l) in pkeys:
                    bg_dup += 1
                    continue
                groups.setdefault((l[5:12], l[77:80], night(jd80(l), l[77:80])), []).append(l)
        seen = set()
        for g, ls in sorted(groups.items()):
            t = fresh()
            for l in ls:
                if key(l) in seen:                                  # the same detection twice in the ITF
                    continue
                seen.add(key(l))
                bg.append((norm(l, t), -1, t))
    rows = planted + bg
    rng.shuffle(rows)
    pl = list(planted)
    rng.shuffle(pl)

    def write80(path, rr):
        open(path, "w").write("\n".join(r[0] for r in rr) + "\n")

    write80(os.path.join(out, "blind.obs80"), rows)
    write80(os.path.join(out, "planted.obs80"), pl)
    with open(os.path.join(out, "blind.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["det_id", "obsTime", "jd_utc", "ra_deg", "dec_deg", "mag", "band", "stn", "trkSub"])
        for n, (l, _, t) in enumerate(rows, 1):
            jd = jd80(l); ra, de = radec(l)
            ts = dt.datetime(2000, 1, 1, 12) + dt.timedelta(days=jd - 2451545.0)
            w.writerow(["blind.obs80:%d" % n, ts.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z", "%.7f" % jd, "%.7f" % ra,
                        "%.7f" % de, l[65:70].strip(), l[70].strip(), l[77:80], t])

    # ---- the null copy: each tracklet moved by its own whole number of days
    by_trk = {}
    for l, oi, t in rows:
        by_trk.setdefault(t, []).append((l, oi))
    per_obj = {}
    for t in sorted(by_trk):
        oi = by_trk[t][0][1]
        if oi >= 0:
            per_obj.setdefault(oi, []).append(t)
    off = {}
    for oi, ts in sorted(per_obj.items()):
        pool = []
        for t in ts:
            if not pool:
                pool = list(NULL_OFFSETS); rng.shuffle(pool)
            jds = [jd80(l) for l, _ in by_trk[t]]
            for k in range(len(pool)):                   # the first offset left that keeps it inside the window
                if t0 <= min(jds) + pool[k] and max(jds) + pool[k] <= t_end:
                    off[t] = pool.pop(k)
                    break
            else:
                off[t] = None
    for t in sorted(by_trk):
        if by_trk[t][0][1] < 0:
            jds = [jd80(l) for l, _ in by_trk[t]]
            k, sg = rng.randint(3, 15), rng.choice((1, -1))
            off[t] = next((s * k for s in (sg, -sg) if t0 <= min(jds) + s * k and max(jds) + s * k <= t_end), None)
    null, nkeys, dropped = [], {}, 0
    for t in sorted(by_trk):
        if off[t] is None:
            dropped += 1
            continue
        for l, oi in by_trk[t]:
            ndec = len(l[23:32].split(".")[1].rstrip()) if "." in l[23:32] else 0
            ndec = max(ndec, 5)
            nl = l[:15] + date80(jd80(l) + off[t], ndec)[:17].ljust(17) + l[32:80]
            null.append((nl, oi, t))
            nkeys[key(nl)] = [oi, t, off[t]]
    rng.shuffle(null)
    write80(os.path.join(out, "null.obs80"), null)

    truth = {"t_end": t_end, "window": window, "seed": seed, "objects": objs,
             "planted_keys": {key(l): oi for l, oi, _ in planted},
             "tracklet_object": {t: oi for l, oi, t in rows if oi >= 0},
             "null_keys": nkeys}
    json.dump(truth, open(os.path.join(out, "truth.json"), "w"))
    fo = [o for o in objs if o["findable"]]
    classes = sorted({o["class"] for o in objs})
    man = {"utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "t_end": t_end, "window": window, "seed": seed,
           "sources": {"sample": [sample, sha(sample)], "background": [bgf, sha(bgf) if bgf else None]},
           "counts": {"planted_rows": len(planted), "background_rows_in_window": bg_in, "background_rows_kept": len(bg),
                      "background_rows_equal_to_planted": bg_dup, "blind_rows": len(rows), "null_rows": len(null),
                      "null_tracklets_left_out": dropped, "tracklets_planted": len({t for _, _, t in planted}),
                      "tracklets_background": len({t for _, _, t in bg}),
                      "objects_with_rows": sum(1 for o in objs if o["rows"]),
                      "findable": {c: sum(1 for o in fo if o["class"] == c) for c in classes}},
           "outputs": {f: sha(os.path.join(out, f)) for f in ("blind.obs80", "blind.csv", "planted.obs80", "null.obs80", "truth.json")}}
    man["counts"]["findable"]["all"] = len(fo)
    json.dump(man, open(os.path.join(out, "manifest.json"), "w"), indent=1)
    print(json.dumps(man["counts"]))


if __name__ == "__main__":
    main()
