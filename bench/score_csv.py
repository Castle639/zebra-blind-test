"""Score any linker's links on a survey-alert copy (rubin_build.py), against the survey's own known-object association.

    python score_csv.py TRUTH.json INPUT.csv LINKS.csv OUT.json [--links-json LINKS.json] [--null] [--designations D.json]
                        [--horizons CACHE.json] [--tol 2.0]

LINKS.csv: link_id, det_id (the input's id column). Each detection's truth: the survey's object id, or 0 (associated with
no known object). A link is
  pure      all its rows one known object's, and every unassociated row in it within --tol arcsec of that object's JPL
            Horizons place (needs --designations: object id -> designation)
  mixed     rows of two or more known objects, or an unassociated row off the object's place
  unverified  one known object's rows and an unassociated row Horizons could not judge (no designation, no ephemeris, or
            an ephemeris that misses the object's own rows in the link by more than --tol)
  unknown   no row of a known object (nothing to score it against)
A known object is findable with 2+ rows on each of 3+ local nights; found with a pure link holding 3+ of its nights.
--null (the copy whose nights' dates were shuffled): a link is a copy when all its rows were moved by one offset (one
night's rows, or nights whose spacing the shuffle kept), else false.
"""
import csv, json, math, os, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score import horizons, sep  # noqa: E402

LON = {"X05": 289.25, "I41": 243.14}     # the stations' longitudes east (ObsCodes), for the local night
STN = ["X05"]


def night(mjd):
    return int(math.floor(mjd + 2400000.5 + LON[STN[0]] / 360.0))


def main():
    a = sys.argv
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d  # noqa: E731
    truth, inp, lcsv, out = json.load(open(a[1])), a[2], a[3], a[4]
    is_null = "--null" in a
    tol = float(opt("--tol", "2.0"))
    des = json.load(open(opt("--designations", None))) if opt("--designations", None) else {}
    hf = opt("--horizons", None)
    cache = json.load(open(hf)) if hf and os.path.exists(hf) else {}
    lj = {x["link_id"]: x for x in json.load(open(opt("--links-json", None)))} if opt("--links-json", None) else {}
    rows = {}
    tai = False
    for r in csv.DictReader(open(inp)):
        tai = "mjd_tai" in r
        STN[0] = r["stn"]
        rows[r["id"]] = (float(r["mjd_tai"] if tai else r["mjd_utc"]), float(r["ra"]), float(r["dec"]))
    dt_utc = 37.0 / 86400 if tai else 0.0
    det = truth["det"]
    noff = {int(k): v for k, v in truth["null_offset"].items()}
    inv = {}                                            # null: the original night of a moved row
    if is_null:
        for n, o in noff.items():
            inv[n + o] = o
    links = defaultdict(list)
    for r in csv.DictReader(open(lcsv)):
        links[r["link_id"]].append(r["det_id"])
    per, found = [], defaultdict(list)
    for lid, ds in links.items():
        L = lj.get(lid, {})
        info = {"link": lid, "rows": len(ds), "nights": len({night(rows[d][0]) for d in ds if d in rows})}
        for f in ("e_fit", "q_fit_au", "v_inf_kms", "rms_arcsec", "n_stations", "provisional"):
            if f in L:
                info[f] = L[f]
        if is_null:
            offs = {inv.get(night(rows[d][0])) for d in ds}
            objs = {det.get(d, 0) for d in ds}
            info.update({"offsets": sorted(o for o in offs if o is not None), "objects": sorted(map(str, objs)),
                         "verdict": "copy" if len(offs) == 1 else "false"})
            per.append(info)
            continue
        lab = defaultdict(list)
        for d in ds:
            lab[det.get(d, 0)].append(d)
        known = [o for o in lab if o]
        if not known:
            info["verdict"] = "unknown"
            per.append(info)
            continue
        verdict, offrows = "mixed", {}
        if len(known) == 1:
            o = known[0]
            verdict = "pure"
            # the arbiter checked first (22:38 UTC 28 Oct, after R05's first reading: 2016 PL293's ephemeris sat 110-130" off
            # its own Rubin rows): Horizons judges an unassociated row only if it places the object's own rows in the link
            # within --tol; otherwise the link is unverified, never mixed. --no-arbiter-check gives the first reading.
            if lab.get(0) and des.get(str(o)) and "--no-arbiter-check" not in a:
                own = []
                for d in lab[o][:2]:
                    m = rows[d]
                    eph = horizons({"desig": des[str(o)], "class": "asteroid"}, STN[0], [m[0] + 2400000.5 - dt_utc], cache)
                    if eph and not isinstance(eph, dict):
                        own.append(sep(m[1], m[2], *eph[0]))
                info["arbiter_own_rows_arcsec"] = [round(x, 2) for x in own]
                if not own or max(own) > tol:
                    verdict = "unverified"
                    offrows = {"arbiter": "Horizons places the object's own rows %s arcsec off" % [round(x, 1) for x in own]}
            for d in (lab.get(0, []) if verdict == "pure" else []):
                if not des.get(str(o)):
                    offrows[d] = "no designation"
                    verdict = "unverified" if verdict == "pure" else verdict
                    continue
                m = rows[d]
                eph = horizons({"desig": des[str(o)], "class": "asteroid"}, STN[0], [m[0] + 2400000.5 - dt_utc], cache)
                if not eph or isinstance(eph, dict):
                    offrows[d] = "no ephemeris"
                    verdict = "unverified" if verdict == "pure" else verdict
                    continue
                s = sep(m[1], m[2], *eph[0])
                offrows[d] = round(s, 2)
                if s > tol:
                    verdict = "mixed"
                    break
            if verdict == "pure":
                nts = defaultdict(int)
                for d in ds:
                    nts[night(rows[d][0])] += 1
                found[o].append(sum(1 for v in nts.values() if v >= 2))
        info.update({"objects": sorted(map(str, known)), "unassociated_rows": len(lab.get(0, [])), "unassociated_check": offrows,
                     "verdict": verdict})
        per.append(info)
    if hf:
        json.dump(cache, open(hf, "w"))
    res = {"input": os.path.basename(inp), "links": len(links)}
    if is_null:
        fl = [x for x in per if x["verdict"] == "false"]
        res.update({"copies": len(per) - len(fl), "false": len(fl),
                    "false_by_nights": {str(n): sum(1 for x in fl if x["nights"] == n) for n in sorted({x["nights"] for x in fl})},
                    "false_provisional": sum(1 for x in fl if x.get("provisional")),
                    "false_on_known_rows_only": sum(1 for x in fl if "0" not in x["objects"])})
    else:
        fo = [o for o, v in truth["objects"].items() if v["nights_2plus"] >= 3]
        fd = [o for o in fo if any(n >= 3 for n in found.get(int(o), []))]
        sc = [x for x in per if x["verdict"] in ("pure", "mixed")]
        res.update({"links_with_known": len(sc), "pure": sum(1 for x in sc if x["verdict"] == "pure"),
                    "mixed": sum(1 for x in sc if x["verdict"] == "mixed"),
                    "purity": round(sum(1 for x in sc if x["verdict"] == "pure") / len(sc), 4) if sc else None,
                    "unverified": sum(1 for x in per if x["verdict"] == "unverified"),
                    "unknown_links": sum(1 for x in per if x["verdict"] == "unknown"),
                    "findable": len(fo), "found": len(fd), "completeness": round(len(fd) / len(fo), 4) if fo else None,
                    "by_nights": {k: {"findable": sum(1 for o in fo if f(truth["objects"][o]["nights_2plus"])),
                                      "found": sum(1 for o in fd if f(truth["objects"][o]["nights_2plus"]))}
                                  for k, f in (("3-4 nights", lambda n: n <= 4), ("5+ nights", lambda n: n >= 5))},
                    "split_objects": sum(1 for v in found.values() if sum(1 for n in v if n >= 3) > 1)})
    json.dump(dict(res, per_link=per), open(out, "w"), indent=1, default=float)
    print(json.dumps(res, default=float))


if __name__ == "__main__":
    main()
