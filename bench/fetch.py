"""Fetch each drawn object's observations from the MPC (get-obs, 80-column), one file an object, resumable.

    python fetch.py OUT/sample.json OUT/obs/

Writes OUT/obs/NNNN.obs80 (NNNN = the object's idx) and OUT/obs/provenance.jsonl (the name asked, the time, rows, sha256).
The MPC is only read; one request at a time.
"""
import hashlib, json, os, re, sys, time, urllib.request

URL = "https://data.minorplanetcenter.net/api/get-obs"


def query_name(o):
    d = o["desig"]
    if o["class"] == "comet":
        d = d.split(" (")[0].strip()
        m = re.match(r"^(\d+P(-[A-Z]+)?)/", d)          # '1P/Halley' -> '1P'; '240P-B/NEAT' -> '240P-B'
        return m.group(1) if m else d
    if d.startswith("("):                               # '(433) Eros' -> '433'
        return d[1:d.index(")")]
    return d


def get_obs(name, tries=5):
    body = json.dumps({"desigs": [name], "output_format": ["OBS80"]}).encode()
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(URL, data=body, method="GET", headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                j = json.loads(r.read())
            return j[0].get("OBS80") or ""
        except Exception as ex:  # noqa: BLE001
            err = ex
            time.sleep(2 ** (k + 1))
    raise RuntimeError("%s: %s" % (name, err))


def main():
    sample, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    S = json.load(open(sample))["objects"]
    prov = open(os.path.join(outdir, "provenance.jsonl"), "a")
    t0, bad = time.time(), 0
    for n, o in enumerate(S):
        f = os.path.join(outdir, "%04d.obs80" % o["idx"])
        if os.path.exists(f):
            continue
        q = query_name(o)
        try:
            txt = get_obs(q)
        except RuntimeError as ex:
            bad += 1
            print("FAIL", o["idx"], q, str(ex)[:200], flush=True)
            continue
        open(f + ".part", "w").write(txt)
        os.replace(f + ".part", f)
        prov.write(json.dumps({"idx": o["idx"], "asked": q, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                               "rows": txt.count("\n"), "sha256": hashlib.sha256(txt.encode()).hexdigest()}) + "\n")
        prov.flush()
        if n % 50 == 0:
            print("fetched %d/%d (%.0f s)" % (n + 1, len(S), time.time() - t0), flush=True)
    print("fetch done: %d objects, %d failed, %.0f s" % (len(S), bad, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
