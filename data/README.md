# data/ — base inputs (vendored)

Read-only base data the pipeline consumes. Vendored here so the repo is
self-contained (nothing is read from the external DeepPBS trees anymore).

| Path | In git? | What |
|---|---|---|
| `folds/` | ✅ yes | Baseline train/valid fold splits + `id.txt` (the benchmark test set). Small text files. |
| `assembly2024/` | ❌ gitignored | The base DeepPBS training feature set (`*.npz`) — large, so kept out of git and copied in once per machine. |

## Populating `assembly2024/` (one-time, on the cluster)

```bash
cp -r /project2/rohs_102/shewchuk/DeepPBS_data/deeppbsmar24/data/assembly2024 \
      "$(git rev-parse --show-toplevel)/data/assembly2024"
```

Override either location via env if you keep them elsewhere:
`ORIG_FOLDS_DIR`, `ORIG_ASSEMBLY_DIR` (see `lib/common.sh`).
