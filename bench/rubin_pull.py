"""Pull a frozen copy of Rubin's public alerts in a sky box and a window, from ALeRCE's TAP (no account).

    python rubin_pull.py MJD0 MJD1 RA0 RA1 DEC0 DEC1 OUTDIR

Three queries, each saved with its provenance (query, time, rows, sha256) in OUTDIR/provenance.jsonl:
  sso.csv   the detections Rubin itself associated with a known solar-system object (detection.sid = 2; their oid is
            Rubin's ssObjectId: the truth)
  dia.csv   every other detection in the box and window (sid = 1; oid = its DIA object)
  diaobj.csv  the DIA objects in the box (oid, ndiasources), so a linker's input can keep the moving-object-like rows
            (a DIA object of 4 sources or fewer)
Resumable (a file present is not asked again); retries with back-off. ALeRCE is only read.
"""
import hashlib, json, os, sys, time, urllib.parse, urllib.request

URL = "https://tap.alerce.online/tap/sync"


def ask(q, fn, prov):
    if os.path.exists(fn):
        return
    for k in range(7):
        try:
            t0 = time.time()
            body = urllib.request.urlopen(URL, data=urllib.parse.urlencode(
                {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": "20000000", "QUERY": q}).encode(),
                timeout=3000).read()
            if body.count(b"\n") < 1 or b"ERROR" in body[:300].upper():
                raise RuntimeError(body[:300])
            break
        except Exception as e:  # noqa: BLE001
            print("retry", os.path.basename(fn), k, str(e)[:200], flush=True)
            time.sleep(15 * 2 ** k)
    else:
        sys.exit("failed: %s" % fn)
    open(fn + ".part", "wb").write(body)
    os.replace(fn + ".part", fn)
    rec = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "url": URL, "query": q, "file": os.path.basename(fn),
           "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(), "rows": body.count(b"\n") - 1,
           "seconds": round(time.time() - t0, 1)}
    prov.write(json.dumps(rec) + "\n")
    prov.flush()
    print(rec["utc"], rec["file"], rec["rows"], rec["seconds"], flush=True)


def main():
    m0, m1, r0, r1, d0, d1 = (float(x) for x in sys.argv[1:7])
    out = sys.argv[7]
    os.makedirs(out, exist_ok=True)
    prov = open(os.path.join(out, "provenance.jsonl"), "a")
    box = "g.ra >= %s AND g.ra < %s AND g.dec >= %s AND g.dec < %s" % (r0, r1, d0, d1)
    win = "g.mjd >= %.6f AND g.mjd < %.6f" % (m0, m1)
    ask("SELECT g.measurement_id, g.oid, g.mjd, g.ra, g.dec, g.band FROM alerce_tap.detection g WHERE g.sid=2 AND %s AND %s"
        % (win, box), os.path.join(out, "sso.csv"), prov)
    ask("SELECT g.measurement_id, g.oid, g.mjd, g.ra, g.dec, g.band FROM alerce_tap.detection g WHERE g.sid=1 AND %s AND %s"
        % (win, box), os.path.join(out, "dia.csv"), prov)
    ask("SELECT g.oid, g.ndiasources, g.firstdiasourcemjdtai, g.lastdiasourcemjdtai FROM alerce_tap.lsst_dia_object g WHERE "
        "g.ra >= %s AND g.ra < %s AND g.dec >= %s AND g.dec < %s" % (r0 - 0.05, r1 + 0.05, d0 - 0.05, d1 + 0.05),
        os.path.join(out, "diaobj.csv"), prov)
    print("done", flush=True)


if __name__ == "__main__":
    main()
