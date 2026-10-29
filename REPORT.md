# Known-object blind test of the Zebra linker: results

29 October 2026

## Summary
| Data set | Detections | Findable objects | Found | Links with a known object: pure / mixed | Purity | False links, null copy |
|---|---|---|---|---|---|---|
| Rubin public alerts (S2) | 78,646 | 842 | **803 (95.4 %)** | 805 / 5 | **99.4 %** | 7 |
| ZTF public alerts (S3) | 44,914 | 1,929 | **1,917 (99.4 %)** | 1,917 / 2 | **99.9 %** | 14 |
| MPC Isolated Tracklet File with planted known objects (S1) | 79,917 | 455 | **416 (91.4 %)** | 440 / 0 | **100 %** | 18 |
| 1I, 2I, 3I at 14 cutoffs, each among the MPC's detections of unnumbered objects of the same weeks (S1-ISO) | 49,404-67,405 per copy | 3 | **3 of 3**, each flagged as interstellar | 0 links mixing an interstellar object with another | crowd 99.1-99.9 % | 0-2 per copy |

All figures are read from the linker's main output, `links.csv`, with the scorer described in `BENCHMARK.md`.

## 1. Test design
- **Truth:** the identifications made by the data providers themselves: the MPC's designations (S1, S1-ISO), Rubin's
  `ssObjectId` (S2), ZTF's `ssnamenr` within 5" (S3). The linker never sees them.
- **Blinding:** every row, known object or background, is rewritten the same way: object numbers, discovery marks, notes and
  publication references blanked; a fresh random tracklet ID per object-station-night that matches no name ever used in the
  ITF; the rows shuffled. CSV rows carry random IDs only. The truth is a separate file.
- **Frozen copies:** every input is stored with its sha256, the sources' fetch times and the queries that produced it
  (`data/*/manifest.json`).
- **Definitions:** an object is *findable* with tracklets of 2+ detections on 3+ local nights in the input. A link is *pure*
  when all its detections belong to one object; a detection without truth counts as that object's only if JPL Horizons places
  it within 2" of the object at every epoch. A link holding a known object's detection that is not pure is *mixed*. An object
  is *found* with a pure link holding 3+ of its nights. Purity = pure / (pure + mixed).
- **Null copies:** the same detections shifted in time so that no real object remains (S1, S1-ISO: each tracklet by its own
  whole number of days, 3-15; S2, S3: the nights' dates permuted, each night's detections moved together). A link whose
  tracklets were all moved by one offset is a *copy* of a real motion; any other link is *false*.

## 2. Software and run conditions
- **Zebra**, one version: sha256 `cb803dca163c82b597562d1b40e2c985b6844625fd5d04003154c60ba7ede599` (the hash of its code
  files' sha256 list, `results/zebra_version.sha256`). Default settings; one command per input; nothing tuned on these copies.
- **Machine:** 4 vCPU (Intel Xeon, 2.80 GHz), 15 GB RAM, Python 3.13; each run alone on the machine, 4 processes.
- **Output:** `links.csv` (one row per detection: `link_id, det_id`) and `links.json` (each link's fitted orbit). Links of
  exactly 3 nights whose free orbit fit is hyperbolic or has chi2/dof above 1.0 (at 0.6" per detection) are written to a
  separate file, `provisional.csv`, pending a fourth night. They are scored the same way and reported separately below; they
  are not counted in the summary.

## 3. Results

### 3.1 Rubin public alerts (S2)
Input: ALeRCE's Rubin alert stream, MJD 61206-61236, RA 330-337.5°, Dec -20 to -15° (6 nights with visits): 13,239 detections
that Rubin associated with 2,744 known objects, and 65,407 detections of DIA objects with 4 or fewer sources (no
association). 78,646 detections, no tracklet IDs (the linker pairs single detections itself).

| | Blind copy | Null copy |
|---|---|---|
| Links | 979 | 7 |
| Pure / mixed / unverified (a) | 805 / 5 / 2 | - |
| Links of unassociated detections only (no truth) | 167 | - |
| Found / findable | 803 / 842 (95.4 %); 3-4 nights 728 / 767, 5+ nights 75 / 75 | - |
| False links | - | 7 (all 3 nights) |
| `provisional.csv` | 58 links: 0 pure, 45 mixed, 1 unverified, 12 without truth | 59 links, all false |
| Run time | 2,304 s | 1,626 s |

(a) One object's detections and an unassociated detection that JPL Horizons could not judge.
Each of the 5 mixed links holds one known object's detections on two nights and unassociated detections on a third, 450-870"
from the object's predicted position. No object is split over two links.

### 3.2 ZTF public alerts (S3)
Input: ALeRCE's ZTF alert stream, MJD 61311.5-61339.5, RA 25-35°, Dec +5 to +15° (near opposition), 17 nights: 44,914
positive-difference detections, 44,847 of them matched by ZTF to a known object within 5".

| | Blind copy | Null copy |
|---|---|---|
| Links | 1,919 | 377 (363 copies) |
| Pure / mixed | 1,917 / 2 | - |
| Found / findable | 1,917 / 1,929 (99.4 %); 3-4 nights 681 / 691, 5+ nights 1,236 / 1,238 | - |
| False links | - | 14 (13 on 3 nights, 1 on 4) |
| `provisional.csv` | 0 links | 111 links: 33 copies, 78 false |
| Run time | 4,559 s | 1,123 s |

Each mixed link joins two known objects (one over 12 nights, one over 4). Each false link in the null copy joins two real
objects' tracklets that were moved by different offsets. No object is split.

### 3.3 MPC Isolated Tracklet File with planted known objects (S1)
Input: the ITF as retrieved 20:12:58 UTC 28 Oct 2026, its 56,647 detections of the 30 days to JD 2461342.342338, with 23,270
of the MPC's own detections of 696 known objects of the same days planted in it (selection seed 65: NEOs, comets near
perihelion, distant objects, main-belt, Hilda, Trojan and other asteroids). 79,917 detections, 23,907 tracklets.

| | Blind copy (B) | Planted detections alone (A) |
|---|---|---|
| Links | 1,367 (440 with a planted detection, 927 of ITF detections only) | 441 |
| Pure / mixed | **440 / 0** | 441 / 0 |
| Found / findable | **416 / 455 (91.4 %)** | 417 / 455 |
| NEOs | 248 / 264 (93.9 %) | 249 / 264 |
| Comets | 44 / 53 (83.0 %) | 44 / 53 |
| Distant (a > 5.5 au) | 50 / 59 (84.7 %) | 50 / 59 |
| Main belt / Hilda / Trojan / other | 35 / 39, 12 / 12, 9 / 9, 18 / 19 | |
| 5+ nights / 3-4 nights | 296 / 301 (98.3 %) / 120 / 154 (77.9 %) | |
| Objects split over 2+ links | 13 | 13 |
| `provisional.csv` | 26 links: 12 with a planted detection, 12 pure, 0 mixed | 11 links, 11 pure |
| Run time | 4,412 s | 3,522 s |

- Every object found in A is found in B except one, which B holds in `provisional.csv`; the 56,647 ITF detections cost no
  object. Two links in B include 4 ITF tracklets, each of the same object by JPL Horizons.
- `provisional.csv` of B holds pure links of 8 objects not found in `links.csv` (6 NEOs, 2 distant objects): 424 of 455
  (93.2 %) are found in the two files together.
- **Orbits** (pure links against the MPC's elements): 5+ nights (292 links) median |Δq|/q 0.022 %, median |Δe| 0.0004;
  3-4 nights (148 links) 0.11 % and 0.0017.
- **Interstellar flag** (e >= 0.95, or a bound fit with chi2/dof > 4; 3+ nights, 2+ stations): raised on 3 of 419 pure links
  of bound objects (0.7 %).
- **Repeatability:** A run a second time gives a byte-identical `links.csv` and `provisional.csv`.
- **Null copy (C):** 105 links: 87 copies, **18 false** (9 on 3 nights, 8 on 4, 1 on 5). 15 of the 18 join one planted
  object's own tracklets at shifted dates, which keep its sky path and rate; 3 join unrelated ITF tracklets. `provisional.csv`:
  100 links, 29 copies, 71 false. Run time 2,215 s.
- Counting nights by calendar date instead of local night: findable 451, found 413 (91.6 %).

### 3.4 Interstellar objects 1I, 2I, 3I (S1-ISO)
Each object's MPC detections (all stations, including precoveries reported after discovery) are cut at its k-th night with a
tracklet; the window is the 30 days to that night. Each cut is linked alone, among the MPC's detections of unnumbered objects
of the same weeks (truth on every detection), and as a null copy of that crowd. "Flagged": the interstellar flag of 3.3.

| Object, k (cutoff, UTC) | Nights (stations) | `links.csv`: alone / among the crowd | `provisional.csv`: alone / among the crowd | Crowd links pure / mixed | Null: false links (`provisional.csv`) |
|---|---|---|---|---|---|
| 1I, 3 (20 Oct 2017 04:27) | 3 (6) | - / - | e 1.23, flagged / e 1.23, flagged | 1,553 / 6 | 2 (2) |
| 1I, 4 (21 Oct 2017 09:23) | 4 (8) | **e 1.17, flagged / e 1.17, flagged** | - | 1,574 / 7 | 0 (5) |
| 1I, 5 (22 Oct 2017 11:38) | 5 (11) | e 1.16, flagged / e 1.16, flagged | - | 1,514 / 7 | 2 (3) |
| 1I, 7 (24 Oct 2017 06:13) | 7 (13) | e 1.18, flagged / e 1.18, flagged | - | 1,638 / 5 | 1 (5) |
| 2I, 6 (1 Sep 2019 03:17) | 3 (3) | - / - | e 2.69, flagged / e 0.66, not flagged (b) | 836 / 4 | 2 (3) |
| 2I, 7 (2 Sep 2019 02:48) | 4 (6) | **e 3.46, flagged / e 3.50, flagged** | - | 932 / 4 | 0 (1) |
| 2I, 8 (3 Sep 2019 01:33) | 5 (7) | e 30.05, flagged / e 30.05, flagged (c) | - | 974 / 3 | 1 (3) |
| 2I, 10 (4 Sep 2019 10:30) | 7 (10) | e 11.37, flagged / e 11.37, flagged | - | 1,006 / 3 | 0 (2) |
| 3I, 3 (8 Jun 2025 09:50) | 3 (1) | - / - | e 6.01, single station / - | 1,422 / 2 | 1 (12) |
| 3I, 4 (18 Jun 2025 10:43) | 4 (1) | - / - | e 6.19, single station / e 6.19, single station | 1,345 / 4 | 1 (9) |
| 3I, 5 (21 Jun 2025 15:15) | 4 (2) | - / - | e 6.19, single station / e 6.19, single station | 1,285 / 4 | 2 (8) |
| 3I, 6 (25 Jun 2025 02:30) | 5 (4) | **e 6.41, flagged / e 6.41, flagged** | - | 1,161 / 8 | 0 (9) |
| 3I, 9 (29 Jun 2025 07:48) | 8 (5) | e 6.11, flagged / e 6.11, flagged | - | 805 / 1 | 0 (6) |
| 3I, 10 (2 Jul 2025 00:18) | 10 (15) | e 6.23, flagged / e 6.23, flagged | - | 776 / 7 | 1 (8) |

(b) The same 17 detections fit a hyperbolic orbit alone and a bound one among the crowd: on a 2.1-day arc the fit has two
acceptable minima. (c) A 5-day arc; e is not constrained.
- In `links.csv`, each object is first linked at its first 4-night link, at the same cutoff alone and among the crowd, and
  flagged from then on: 1I two days after its discovery, 2I three days after, 3I on 25 June 2025, six days before its
  discovery of 1 July (from precoveries). The MPC's eccentricities: 1I 1.20, 2I 3.36, 3I 6.1.
- In `provisional.csv`, each is linked at its first 3-night cutoff: 1I the day after its discovery; 3I from ZTF (station I41)
  detections to 8 June 2025, 23 days before its discovery.
- No link mixes an interstellar object with another object. The crowd's own links are 99.1-99.9 % pure in every copy.
- Run times: alone 9-32 s; among the crowd 423-1,004 s; null copies 263-628 s.

## 4. Files
- `BENCHMARK.md`: the copies, what is hidden, the definitions, how to score any linker, how to build a new copy.
- `data/s1`, `data/s2`, `data/s3`, `data/iso`: the inputs, null copies, truth and manifests; `data/ObsCodes.html`, the MPC's
  observatory codes; the JPL Horizons positions used by the scorer.
- `bench/`: the copy builders and the scorer (`score.py`, `score_csv.py`, `test_score.py`; standard-library Python).
- `results/`: for every input, Zebra's `links.csv`, `links.json`, `provisional.csv`, `provisional.json`, the scores
  (`*score.json`) and the run summary (input sha256, linker file hashes, settings, run time).
