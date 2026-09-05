# Exploratory Milo analyses

Run from the repository root after configuring the external data directory.
These analyses are not part of the manuscript pipeline.

- `miloR-da.R` builds neighborhoods in PCA space.
- `module-score-milo-da.R` builds neighborhoods from four standardized marker
  module scores. Its default primary settings are `k = 60` and `prop = 0.04`.

```sh
Rscript exploratory/milo/miloR-da.R
Rscript exploratory/milo/module-score-milo-da.R
```

Both use Mouse × Condition samples for differential-abundance inference and
retain their existing external output paths under `degs/mg_selected/`.
Set `ESPI_OVERWRITE=true` only when deliberately replacing results. The module
score script exposes `ESPI_MODULE_SCORE_MILO_K`, `ESPI_MODULE_SCORE_MILO_PROP`,
and `ESPI_MODULE_SCORE_MILO_OUTPUT_DIR` for its existing run settings.

The [neurogenic-proportion sensitivity analysis](../neurogenic-proportion/notes.md)
consumes the primary module-score table. Moving these scripts does not change
that dependency or the saved results.
