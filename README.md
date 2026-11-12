# SynthNetGen_Mod

## Setup
```bash
chmod +x build.sh
./build.sh
```

## Usage

Generate a starting graph automatically:

```bash
./rpll.sh <k:4|5> <target> <output> <RSEAS_rounds> <RMSE_rounds> <r>
```

Or supply an existing starting graph:

```bash
./rpll.sh <k:4|5> <target> <output> <synth_base> <RSEAS_rounds> <RMSE_rounds> <r>
```

Use `./rpllb.sh` in either command to estimate graphlet frequencies with BLANT instead of ORCA.
It still uses ORCA for per-node orbit counts. Use `./rplld.sh` to combine BLANT
frequency estimates with sampled binary node-orbit presence from `sampled_orbits`,
without requiring ORCA at runtime. All three scripts accept the same arguments.

- `k`: graphlet size, either 4 or 5.
- `target`: target graph edge-list file.
- `output`: destination for the generated graph.
- `synth_base`: optional starting graph with the same node count as the cleaned target.
- `RSEAS_rounds`: number of RSEAS optimization rounds; may be zero.
- `RMSE_rounds`: number of subsequent RMSE optimization rounds; may be zero.
- `r`: baseline results are computed every `r` rounds.

For parallel workers, use

```bash
PARALLEL_JOBS=4 ./rpllb.sh 4 target.el output.el 12 3 3
```

`rplld.sh` defaults to one worker and 10,000,000 orbit samples per node-table
preparation. Set `ORBIT_SAMPLES` to change the sample budget and optionally set
`ORBIT_SEED` for reproducible orbit sampling:

```bash
ORBIT_SAMPLES=1000000 ORBIT_SEED=42 ./rplld.sh 4 target.el output.el 12 3 3
```

Set `PARALLEL_JOBS` and `ROUND_JOBS` to change its probe and round-worker limits.
