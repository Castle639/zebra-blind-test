# Known-object blind test for linkers

- `REPORT.md`: the results of the Zebra linker on this test.
- `BENCHMARK.md`: the test design, the definitions, how to score any linker's links, and how to build a new copy.
- `data/`: the frozen inputs (the MPC Isolated Tracklet File with planted known objects; 1I, 2I and 3I at 14 cutoffs; Rubin and
  ZTF public alerts), their null copies, the truth files and the manifests.
- `bench/`: the copy builders and the scorer (Python 3.11-3.13, standard library only).
- `results/`: Zebra's links, scores and run summaries for every input.
- `SHA256SUMS`: the sha256 of every file in this folder.
