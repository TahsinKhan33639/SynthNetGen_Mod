# SynthNetGen_Mod

## Setup
```bash
chmod +x build.sh
./build.sh
```

## Usage
```bash
./rpll.sh <k:4|5> <target> <output> <RSEAS_optimizing_rounds> <RMSE_optimizing_round> <r>
```
or
```bash
./rpll.sh <k:4|5> <target> <output> <synth_base> <RSEAS_optimizing_rounds> <RMSE_optimizing_round> <r>
```

Where `k` is the graphlet size and the baseline resuls are computed every `r` rounds.

Use `rpllb.sh` instead of `rpll.sh` to use `blant` instead or `orca`.
