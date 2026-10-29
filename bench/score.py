"""Score any linker's links against a blind copy's truth.

    python score.py TRUTH.json INPUT.obs80 LINKS.csv OUT.json [--links-json LINKS.json] [--null] [--horizons CACHE.json] [--tol 2.0]
                    [--obscodes ObsCodes.html]

LINKS.csv: one row a detection, with columns link_id and det_id; det_id is "<input file name>:<line>" (zlink writes it
so; blind.csv carries it for any other linker). INPUT.obs80 is the file the linker was given (blind.obs80,
planted.obs80 or null.obs80). --links-json (zlink's links.json) adds each link's fitted orbit for the accuracy table and
the interstellar flag (e >= 0.95, or a bound fit with chi2/dof > 4; 3+ nights, 2+ stations).

Definitions:
  tracklet   the rows of one tracklet ID (columns 6-12) at one station on one local night
  findable   an object with tracklets of 2+ rows on 3+ local nights in the window
  pure       a link whose planted rows all belong to one object, and whose background tracklets (if any) each lie within
             --tol arcsec of that object's JPL Horizons place on every row (else the link is mixed)
  mixed      any other link holding a planted row
  found      a findable object with a pure link holding 3+ of its local nights
  split      a found object with 2+ such links
  completeness = found / findable (by class); purity = pure / (pure + mixed)
  --null: every link is a copy (all its tracklets one object's, or all background, moved by one offset: a real motion
  moved in time, not a link by chance) or false (any other); the false links are the answer to "how many links does a
  copy with no real object in it give".
"""
import csv, json, math, os, sys, time, urllib.parse, urllib.request
from collections import defaultdict

E_SEL, BOUND_SEL = 0.95, 4.0     # the interstellar flag: e >= E_SEL, or a bound fit with chi2/dof > BOUND_SEL


def jd80(l):
    y, m, d = int(l[15:19]), int(l[20:22]), float(l[23:32])
    a = (14 - m) // 12; yy = y + 4800 - a; mm = m + 12 * a - 3
    return int(d) + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045 - 0.5 + (d - int(d))


def radec(l):
    h, mi, s = float(l[32:34]), float(l[35:37]), float(l[38:44])
    sg = -1.0 if l[44] == "-" else 1.0
    dd, dm, ds = float(l[45:47]), float(l[48:50]), float(l[51:56])
    return 15.0 * (h + mi / 60 + s / 3600), sg * (dd + dm / 60 + ds / 3600)


def key(l):
    return l[15:56] + l[77:80]


def horizons_cmd(o):
    d = o["desig"]
    if o["class"] == "comet":
        import re
        d = d.split(" (")[0].strip()
        m = re.match(r"^(\d+P(-[A-Z]+)?)/", d)
        d = m.group(1) if m else d
        return "DES=%s;CAP;NOFRAG" % d if "-" not in d else "DES=%s;CAP" % d
    if o["class"] == "interstellar":
        return {"1I": "A/2017 U1", "2I": "C/2019 Q4", "3I": "C/2025 N1"}.get(d, d)
    if d.startswith("("):
        return "%s;" % d[1:d.index(")")]
    if d.isdigit():                       # a numbered asteroid, as ZTF and Rubin name it
        return "%s;" % d
    return "DES=%s;" % d


def horizons(o, stn, jds, cache):
    k = "%s|%s|%s" % (horizons_cmd(o), stn, ",".join("%.6f" % j for j in jds))
    if k in cache:
        return cache[k]
    q = {"format": "json", "COMMAND": "'%s'" % horizons_cmd(o), "EPHEM_TYPE": "OBSERVER", "CENTER": "'%s@399'" % stn,
         "TLIST_TYPE": "JD", "TIME_TYPE": "UT", "TLIST": "'%s'" % " ".join("%.6f" % j for j in jds), "QUANTITIES": "'1'",
         "ANG_FORMAT": "DEG", "EXTRA_PREC": "YES", "CSV_FORMAT": "YES"}
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(q)
    for n in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                res = json.loads(r.read())["result"]
            break
        except Exception:  # noqa: BLE001
            time.sleep(2 ** (n + 1))
    else:
        return None
    if "$$SOE" not in res:
        cache[k] = {"error": res[:300]}
        return cache[k]
    out = []
    for l in res.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines():
        f = [x.strip() for x in l.split(",")]
        nums = [x for x in f[1:] if x not in ("", "*", "C", "N", "A", "m")]
        out.append((float(nums[0]), float(nums[1])))
    cache[k] = out
    return out


def sep(ra1, d1, ra2, d2):
    r = math.radians
    c = math.sin(r(d1)) * math.sin(r(d2)) + math.cos(r(d1)) * math.cos(r(d2)) * math.cos(r(ra1 - ra2))
    return math.degrees(math.acos(max(-1.0, min(1.0, c)))) * 3600.0


def cal(l, _lon={}):
    """A row's calendar night: floor(JD + signed longitude / 360), the same label for one evening at every station."""
    if not _lon:
        _lon.update(obscodes())
    L = _lon.get(l[77:80], 0.0)
    return math.floor(jd80(l) + (((L + 180.0) % 360.0) - 180.0) / 360.0)


OBSCODES = ["data/ObsCodes.html" if os.path.exists("data/ObsCodes.html") else "ObsCodes.html"]   # the MPC's observatory codes


def obscodes():
    lon = {}
    for l in open(OBSCODES[0], errors="replace"):
        if len(l) > 30 and l[3] == " ":
            try:
                lon[l[0:3]] = float(l[4:13])
            except ValueError:
                pass
    return lon


def main():
    a = sys.argv
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d  # noqa: E731
    truth, inp, lcsv, out = json.load(open(a[1])), a[2], a[3], a[4]
    is_null = "--null" in a
    hf = opt("--horizons", None)
    tol = float(opt("--tol", "2.0"))
    cache = json.load(open(hf)) if hf and os.path.exists(hf) else {}
    lj = {x["link_id"]: x for x in json.load(open(opt("--links-json", None)))} if opt("--links-json", None) else {}
    OBSCODES[0] = opt("--obscodes", OBSCODES[0])
    lon = obscodes()
    objs = {o["idx"]: o for o in truth["objects"]}
    lines = [l.rstrip("\n") for l in open(inp)]
    base = os.path.basename(inp)
    tobj = {k: int(v) for k, v in truth["planted_keys"].items()}
    nkeys = truth["null_keys"]
    # findable is read from the very file the linker was given: an object's rows there, as tracklets of 2+ rows on 3+ nights
    tk = defaultdict(int)
    for l in lines:
        if len(l) >= 80 and key(l) in tobj:
            tk[(tobj[key(l)], l[5:12], l[77:80], math.floor(jd80(l) + lon.get(l[77:80], 0.0) / 360.0))] += 1
    onights = defaultdict(set)
    for (oi, t, s, nt), c in tk.items():
        if c >= 2:
            onights[oi].add(nt)
    for o in objs.values():
        o["nights"] = len(onights.get(o["idx"], ()))
        o["findable"] = o["nights"] >= 3
        o["in_input"] = any(k[0] == o["idx"] for k in tk) if not is_null else None
    links = defaultdict(list)
    n_unread = 0
    cal_obj = defaultdict(set)                       # calendar nights (signed longitude) of each object's 2+ row tracklets
    for l in lines:
        if len(l) >= 80 and key(l) in tobj and tk[(tobj[key(l)], l[5:12], l[77:80], math.floor(jd80(l) + lon.get(l[77:80], 0.0) / 360.0))] >= 2:
            cal_obj[tobj[key(l)]].add(cal(l))
    for r in csv.DictReader(open(lcsv)):
        d = r["det_id"]
        if ":" in d and d.rsplit(":", 1)[0] == base:
            links[r["link_id"]].append(lines[int(d.rsplit(":", 1)[1]) - 1])
        else:
            n_unread += 1
    per = []
    found = defaultdict(list)
    cal_found = defaultdict(list)
    for lid, ls in links.items():
        trk = defaultdict(list)
        for l in ls:
            trk[(l[5:12], l[77:80], math.floor(jd80(l) + lon.get(l[77:80], 0.0) / 360.0))].append(l)
        L = lj.get(lid, {})
        info = {"link": lid, "rows": len(ls), "tracklets": len(trk), "nights": len({k[2] for k in trk}),
                "calendar_nights": len({cal(l) for l in ls})}
        for f in ("e_fit", "q_fit_au", "v_inf_kms", "rms_arcsec", "free_fit_chi2_dof", "bound_fit_chi2_dof", "n_stations", "provisional"):
            if f in L:
                info[f] = L[f]
        if is_null:
            tags = {(nkeys[key(v[0])][0] if nkeys[key(v[0])][0] >= 0 else "bg", nkeys[key(v[0])][2]) for v in trk.values()}
            objs_in = {t[0] for t in tags}
            offs = {t[1] for t in tags}
            info.update({"objects": sorted(map(str, objs_in)), "offsets": sorted(offs),
                         "verdict": "copy" if len(offs) == 1 and len(objs_in) == 1 else "false"})
            per.append(info)
            continue
        mine = {k: tobj[key(v[0])] for k, v in trk.items() if key(v[0]) in tobj}
        if not mine:
            info["verdict"] = "background only"
            per.append(info)
            continue
        oset = set(mine.values())
        bgt = [k for k in trk if k not in mine]
        bgv = {}
        if len(oset) == 1 and bgt:
            o = objs[next(iter(oset))]
            for k in bgt:
                rr = trk[k]
                eph = horizons(o, k[1], [jd80(l) for l in rr], cache)
                if not eph or isinstance(eph, dict) or len(eph) != len(rr):
                    bgv["%s|%s|%d" % k] = {"same": False, "why": "no ephemeris"}
                    continue
                s = [sep(*radec(l), *e) for l, e in zip(rr, eph)]
                bgv["%s|%s|%d" % k] = {"same": max(s) <= tol, "max_sep_arcsec": round(max(s), 2)}
        pure = len(oset) == 1 and all(v["same"] for v in bgv.values())
        nts = defaultdict(set)
        for k, oi in mine.items():
            if len(trk[k]) >= 2:
                nts[oi].add(k[2])
        if pure:
            oi = next(iter(oset))
            for k in bgt:
                if len(trk[k]) >= 2:
                    nts[oi].add(k[2])
            found[oi].append(len(nts[oi]))
            cal_found[oi].append(len({cal(l) for k in trk if len(trk[k]) >= 2 for l in trk[k]}))
        info.update({"objects": sorted(oset), "classes": sorted({objs[o]["class"] for o in oset}),
                     "planted_tracklets": len(mine), "background_tracklets": len(bgt), "background_verdict": bgv,
                     "verdict": "pure" if pure else "mixed"})
        if L:
            info["interstellar_flag"] = bool((L.get("e_fit", 0) >= E_SEL or (L.get("bound_fit_chi2_dof") or 0) > BOUND_SEL)
                                          and info["nights"] >= 3 and L.get("n_stations", 0) >= 2)
        per.append(info)
    if hf:
        json.dump(cache, open(hf, "w"))
    res = {"input": base, "links": len(links), "detections_not_in_input": n_unread}
    if is_null:
        res.update({"copies": sum(1 for x in per if x["verdict"] == "copy"),
                    "false": sum(1 for x in per if x["verdict"] == "false"),
                    "false_by_nights": dict(sorted(defaultdict(int, {n: sum(1 for x in per if x["verdict"] == "false" and x["nights"] == n)
                                                                     for n in {x["nights"] for x in per}}).items())),
                    "false_provisional": sum(1 for x in per if x["verdict"] == "false" and x.get("provisional")),
                    "false_hyperbolic": sum(1 for x in per if x["verdict"] == "false" and (x.get("e_fit") or 0) >= 1),
                    "false_with_background": sum(1 for x in per if x["verdict"] == "false" and "bg" in x["objects"])})
    else:
        classes = sorted({o["class"] for o in objs.values()})
        by = {}
        for c in classes + ["all"]:
            fo = [o for o in objs.values() if o["findable"] and (c == "all" or o["class"] == c)]
            fd = [o for o in fo if any(n >= 3 for n in found.get(o["idx"], []))]
            by[c] = {"findable": len(fo), "found": len(fd), "completeness": round(len(fd) / len(fo), 4) if fo else None,
                     "missed": sorted(o["desig"] for o in fo if o not in fd)}
        fo = [o for o in objs.values() if o["findable"]]
        nb = {"3-4 nights": [o for o in fo if o["nights"] <= 4], "5+ nights": [o for o in fo if o["nights"] >= 5]}
        withp = [x for x in per if x["verdict"] in ("pure", "mixed")]
        npure = sum(1 for x in withp if x["verdict"] == "pure")
        acc = []
        for x in withp:
            if x["verdict"] == "pure" and "q_fit_au" in x:
                o = objs[x["objects"][0]]
                if o.get("q") and o.get("e") is not None and o["class"] != "interstellar":
                    acc.append({"desig": o["desig"], "class": o["class"], "nights": x["nights"], "q": o["q"], "e": o["e"],
                                "q_fit": x["q_fit_au"], "e_fit": x["e_fit"], "interstellar_flag": x.get("interstellar_flag")})
        med = lambda v: sorted(v)[len(v) // 2] if v else None  # noqa: E731
        res.update({"links_with_planted": len(withp), "pure": npure, "mixed": len(withp) - npure,
                    "purity": round(npure / len(withp), 4) if withp else None,
                    "background_only_links": sum(1 for x in per if x["verdict"] == "background only"),
                    "links_taking_background_tracklets": sum(1 for x in withp if x["background_tracklets"]),
                    "background_tracklets_same_object": sum(1 for x in withp for v in x["background_verdict"].values() if v["same"]),
                    "completeness": by,
                    "by_nights": {k: {"findable": len(v), "found": sum(1 for o in v if any(n >= 3 for n in found.get(o["idx"], [])))}
                                  for k, v in nb.items()},
                    "split_objects": sum(1 for v in found.values() if sum(1 for n in v if n >= 3) > 1),
                    "calendar_nights": {"links_under_3": sum(1 for x in per if x["calendar_nights"] < 3),
                                        "planted_links_under_3": sum(1 for x in withp if x["calendar_nights"] < 3),
                                        "findable": sum(1 for o in fo if len(cal_obj.get(o["idx"], ())) >= 3),
                                        "found": sum(1 for o in fo if len(cal_obj.get(o["idx"], ())) >= 3
                                                     and any(n >= 3 for n in cal_found.get(o["idx"], [])))},
                    "accuracy": {"5+ nights": {"links": sum(1 for x in acc if x["nights"] >= 5),
                                               "median_dq_over_q": med([abs(x["q_fit"] - x["q"]) / x["q"] for x in acc if x["nights"] >= 5]),
                                               "median_de": med([abs(x["e_fit"] - x["e"]) for x in acc if x["nights"] >= 5])},
                                 "3-4 nights": {"links": sum(1 for x in acc if x["nights"] < 5),
                                                "median_dq_over_q": med([abs(x["q_fit"] - x["q"]) / x["q"] for x in acc if x["nights"] < 5]),
                                                "median_de": med([abs(x["e_fit"] - x["e"]) for x in acc if x["nights"] < 5])}},
                    "interstellar_flag": {"bound_objects_pure_links": sum(1 for x in acc if x["e"] < 0.95),
                              "bound_objects_flagged": sorted({x["desig"] for x in acc if x["e"] < 0.95 and x["interstellar_flag"]}),
                              "interstellar_flagged": sorted({objs[x["objects"][0]]["desig"] for x in withp if x["verdict"] == "pure"
                                                               and objs[x["objects"][0]]["class"] == "interstellar" and x.get("interstellar_flag")})}})
    json.dump(dict(res, per_link=per), open(out, "w"), indent=1, default=float)
    print(json.dumps(res, default=float)[:3000])


if __name__ == "__main__":
    main()
