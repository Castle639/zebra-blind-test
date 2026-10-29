"""Pull a frozen copy of ZTF's public alerts in a sky box and a window, from ALeRCE's TAP (no account), one day a query.

    python ztf_pull.py MJD0 MJD1 RA0 RA1 DEC0 DEC1 OUTDIR [--step DAYS] [--no-obj]

Three queries a day, each saved with its provenance (query, time, rows, sha256) in OUTDIR/provenance.jsonl:
  det_<mjd>.csv   every ZTF alert detection in the box (detection.sid = 0): measurement_id, oid, mjd, ra, dec, band, with
                  ztf_detection's rb, drb, magpsf, isdiffpos
  ss_<mjd>.csv    ZTF's own known-object match of each of those detections (ztf_ss: ssnamenr, ssdistnr): the truth
  obj_<mjd>.csv   each of their ALeRCE objects' ndethist (ztf_object), so a linker's input can keep the moving-object-like rows
Resumable; retries with back-off (rubin_pull.ask). ALeRCE is only read. (The whole window in one joined query was dropped by
the server, 21:13 UTC 28 Oct 2026.)
"""
import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rubin_pull import ask  # noqa: E402


def main():
    m0, m1, r0, r1, d0, d1 = (float(x) for x in sys.argv[1:7])
    out = sys.argv[7]
    step = float(sys.argv[sys.argv.index("--step") + 1]) if "--step" in sys.argv else 1.0
    no_obj = "--no-obj" in sys.argv          # (21:33 UTC 28 Oct: 3-4 min a query; ztf_build then counts each oid's rows in the pull)
    os.makedirs(out, exist_ok=True)
    prov = open(os.path.join(out, "provenance.jsonl"), "a")
    a = m0
    while a < m1 - 1e-9:
        b = min(a + step, m1)
        w = "g.sid=0 AND g.mjd >= %.6f AND g.mjd < %.6f AND g.ra >= %s AND g.ra < %s AND g.dec >= %s AND g.dec < %s" % (a, b, r0, r1, d0, d1)
        ask("SELECT g.measurement_id, g.oid, g.mjd, g.ra, g.dec, g.band, z.rb, z.drb, z.magpsf, z.isdiffpos FROM alerce_tap.detection g "
            "JOIN alerce_tap.ztf_detection z ON z.measurement_id=g.measurement_id WHERE " + w, os.path.join(out, "det_%.1f.csv" % a), prov)
        ask("SELECT s.measurement_id, s.ssnamenr, s.ssdistnr, s.ssmagnr FROM alerce_tap.ztf_ss s JOIN alerce_tap.detection g "
            "ON g.measurement_id=s.measurement_id WHERE " + w, os.path.join(out, "ss_%.1f.csv" % a), prov)
        if not no_obj:
          ask("SELECT DISTINCT o.oid, o.ndethist, o.ncovhist FROM alerce_tap.ztf_object o JOIN alerce_tap.detection g ON g.oid=o.oid WHERE " + w,
              os.path.join(out, "obj_%.1f.csv" % a), prov)
        a = b
    print("done", flush=True)


if __name__ == "__main__":
    main()
