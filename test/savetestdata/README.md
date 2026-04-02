# Save Baseline Test Data

This directory generates deterministic baseline data for top-level libjupitermag APIs.

## What It Saves

- `../data/field_baseline.csv`
  - Internal (`InternalField`), external (`Con2020Field`), and combined (`ModelField`) vectors at fixed points.
- `../data/trace_summary.csv`
  - Per-trace footprint outputs (ionosphere/surface/equatorial footprint positions and lat/lon values).

## Build And Run

```bash
cd test/savetestdata
make
```

This creates/overwrites files in `test/data`.

## Notes

- Uses only public top-level APIs from `jupitermag.h`.
- Intended as baseline generation before API refactor work.
