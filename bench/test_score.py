"""The scorer checked against links whose answer is known before it runs (a control on the ruler itself).

    python test_score.py TRUTH.json INPUT.obs80 LINKS.csv WORKDIR

From a real run's links.csv: (1) unchanged: the scorer's counts; (2) one planted row swapped between two pure links of
different objects: exactly 2 more mixed links, 2 fewer pure; (3) every planted link cut to its first 2 nights: no object
found; (4) one link duplicated under a new id: one more split object (if its object was found once). Prints PASS or FAIL.
"""
import csv, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def score(truth, inp, rows, work, tag):
    p = os.path.join(work, "links_%s.csv" % tag)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["link_id", "det_id"])
        w.writerows(rows)
    out = os.path.join(work, "score_%s.json" % tag)
    subprocess.run([sys.executable, os.path.join(HERE, "score.py"), truth, inp, p, out], check=True, capture_output=True)
    return json.load(open(out))


def main():
    truth, inp, lcsv, work = sys.argv[1:5]
    os.makedirs(work, exist_ok=True)
    rows = [(r["link_id"], r["det_id"]) for r in csv.DictReader(open(lcsv))]
    s0 = score(truth, inp, rows, work, "0")
    pure = [x for x in s0["per_link"] if x["verdict"] == "pure" and x["background_tracklets"] == 0]
    a, b = pure[0], next(x for x in pure if x["objects"] != pure[0]["objects"])
    ra = next(i for i, r in enumerate(rows) if r[0] == a["link"])
    rb = next(i for i, r in enumerate(rows) if r[0] == b["link"])
    sw = list(rows)
    sw[ra], sw[rb] = (rows[ra][0], rows[rb][1]), (rows[rb][0], rows[ra][1])
    s1 = score(truth, inp, sw, work, "swap")
    T = json.load(open(truth))
    lines = [l.rstrip("\n") for l in open(inp)]
    import math
    import score as S
    lon = S.obscodes()
    nt = lambda l: math.floor(S.jd80(l) + lon.get(l[77:80], 0.0) / 360.0)  # noqa: E731  (the scorer's own night label)
    cut, nights = [], {}
    for lid, d in rows:
        nights.setdefault(lid, set()).add(nt(lines[int(d.rsplit(":", 1)[1]) - 1]))
    for lid, d in rows:
        if nt(lines[int(d.rsplit(":", 1)[1]) - 1]) in sorted(nights[lid])[:2]:
            cut.append((lid, d))
    s2 = score(truth, inp, cut, work, "cut")
    once = next(x for x in pure if x["nights"] >= 3)
    dup = rows + [("DUP", d) for lid, d in rows if lid == once["link"]]
    s3 = score(truth, inp, dup, work, "dup")
    ok = [s1["mixed"] == s0["mixed"] + 2 and s1["pure"] == s0["pure"] - 2,
          s2["completeness"]["all"]["found"] == 0,
          s3["split_objects"] >= s0["split_objects"] + 1 and s3["completeness"]["all"]["found"] == s0["completeness"]["all"]["found"]]
    print(json.dumps({"base": [s0["pure"], s0["mixed"], s0["completeness"]["all"]["found"], s0["split_objects"]],
                      "swap": [s1["pure"], s1["mixed"]], "cut_found": s2["completeness"]["all"]["found"],
                      "dup_split": s3["split_objects"], "checks": ok}))
    print("score.py control:", "PASS" if all(ok) else "FAIL")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.path.insert(0, HERE)
    sys.exit(main())
