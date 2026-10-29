# A known-object blind test for linkers: the copies, the scorer, and how to make your own

A benchmark for linking: real detections of known asteroids, comets and interstellar objects, their names taken off, mixed
into real survey data. A linker is given the file and nothing else, and its links are scored against the identifications the
MPC, Rubin and ZTF made themselves. Every input is frozen and hashed here, and the truth is a separate file. The scorer reads
any linker's output (one row per detection: `link_id, det_id`), so the same copy scores any linker alike.

## The copies (`data/`)
| Copy | What a linker is given | Truth (never in the linker's file) |
|---|---|---|
| **S1, the ITF stream** (`data/s1/`) | `blind.obs80` (MPC 80-column; `blind.csv` the same rows as CSV): the MPC's Isolated Tracklet File as fetched 20:12:58 UTC 28 Oct 2026, its 30 days to JD 2461342.342338, with the MPC's own rows of 696 known objects of the same days planted in it; `planted.obs80`, the planted rows alone | the MPC's designations (`truth.json`), each object's MPC rows from get-obs |
| **S1-ISO** (`data/iso/`) | 1I, 2I and 3I, each cut at its k-th night with a tracklet (the window the 30 days to that night): alone (`<ISO>_k<k>_alone.obs80`) and among the MPC's rows of unnumbered objects of the same weeks (`<ISO>_k<k>.obs80`) | the MPC's designations, on every row (`<ISO>_truth.json`) |
| **S2, Rubin's public alerts** (`data/s2/`) | `blind.csv`: every Rubin alert detection of 30 days in a 37.5 deg² box that Rubin tied to a known object, and every other detection of a DIA object with 4 sources or fewer | Rubin's own ssObjectId |
| **S3, ZTF's public alerts** (`data/s3/`) | `blind.csv`: ZTF's alert detections of 28 days in a 98 deg² box, the known-object matches and the moving-object-like rest | ZTF's own ssnamenr (within 5") |

Each copy has a **null copy** (`null.obs80`, `null.csv`, `<ISO>_k<k>_null.obs80`): the same rows moved in time so that no
real object is left in it. S1 and S1-ISO move every tracklet by its own whole number of days, 3 to 15. S2 and S3 shuffle the
nights' dates among themselves, each night's rows moved together. Every link a linker makes on a null copy is counted.

Each copy's `manifest.json` holds the sha256 of every input and output, the sources' fetch times and the queries that made it.
`data/ObsCodes.html` is the MPC's list of observatory codes used for every local night (sha256
`58454e19015f622472a67928d517ca081abfa13c987a211dbf4533e293226fe4`). `data/s1`, `s2` and `s3` each hold the JPL Horizons
positions the scorer used (`horizons_cache.json`), so the scores here can be reproduced offline.

## What is hidden
Every row is rewritten the same way, planted or background. In the 80-column files, columns 1-5 (number), 13 (discovery
mark), 14 (note 1), 57-65 and 73-77 (publication reference) are blank, and a comet's N or T magnitude flag is blank. Columns
6-12 hold a fresh random tracklet ID, one per object-station-night or per original trkSub-station-night. No ID is any name
that ever appeared in the ITF or in the objects' rows (`bench/check_blind.py`). The rows are shuffled. In the CSVs every
row's id is random, with no oid, measurement or candidate id, and no survey name.

## Definitions (`bench/score.py`, `bench/score_csv.py`)
- **tracklet:** a station's rows of one tracklet ID on one local night (S1, S1-ISO), or the linker's own pairing of single
  rows (S2, S3)
- **findable:** an object with tracklets of 2+ rows on 3+ local nights in the file the linker was given
- **pure link:** all its rows one object's. A row with no truth (an ITF or unassociated row) counts as that object's only
  when JPL Horizons puts it within 2" of the object on every row
- **mixed link:** any other link holding a known object's row
- **found:** a findable object with a pure link holding 3+ of its nights; **split:** found by 2+ such links
- **completeness** = found / findable, by class; **purity** = pure / (pure + mixed)
- **null copy:** a link is a **copy** when all its tracklets were moved by one offset (a real motion moved in time), else **false**

## Score a linker's links (Python 3.11-3.13, the standard library only; run from this folder)
Give the linker `blind.obs80` or `blind.csv`, and write its links as `link_id,det_id`: in the 80-column copies det_id is
`<file name>:<line>` (the CSV's `det_id` column carries it), in S2 and S3 it is the CSV's `id`.
```bash
gunzip -k data/*/*.gz
python bench/score.py data/s1/truth.json data/s1/blind.obs80 my/links.csv my/score.json --horizons data/s1/horizons_cache.json
python bench/score.py data/s1/truth.json data/s1/null.obs80 my/null_links.csv my/null_score.json --null
python bench/score_csv.py data/s2/truth.json data/s2/blind.csv my/links.csv my/score.json \
  --designations data/s2/designations.json --horizons data/s2/horizons_cache.json     # S3 the same with s3
python bench/score.py data/iso/1I_truth.json data/iso/1I_k3.obs80 my/links.csv my/score.json
```
`--links-json` (an orbit per link, in the form of `results/*/*links.json`) adds the orbit table and the interstellar flag.
The scorer reads the station list from `data/ObsCodes.html` (`--obscodes` for another path) and asks JPL Horizons only for
positions not in the cache. Zebra's scores in `results/` are reproduced by the same lines on its own `links.csv`.
`bench/test_score.py` checks the scorer itself: a row swapped between two links reads 2 mixed, links cut to 2 nights find
nothing, and a duplicated link reads one more split.

## Make your own copy (your objects, your seed)
```bash
python bench/draw.py MPCORB.DAT.gz CometEls.txt my/sample.json --seed <yours> --since <YYYYMMDD> --counts NEO=100,Distant=all
python bench/fetch.py my/sample.json my/obs                      # the MPC's get-obs
python bench/build.py my/sample.json my/obs my/copy --t-end <JD> --seed <yours> --background itf.txt.gz
python bench/check_blind.py my/copy my/sample.json my/obs --background itf.txt.gz
```
Keep `my/copy/truth.json` to yourself, send us `blind.obs80` alone, and score the links we send back.
