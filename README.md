# TB Screening Cost-Effectiveness Analysis — UK Migrants

Cost-effectiveness analysis of 43 tuberculosis screening strategies for migrants entering the UK. Coupled decision tree–Markov cohort model in R.

## Requirements

- R ≥ 4.4
- Packages: `dplyr`, `ggplot2`, `patchwork`, `readxl`, `scales`, `tibble`, `parallel`

Install all at once:
```r
install.packages(c("dplyr", "ggplot2", "patchwork", "readxl", "scales", "tibble"))
```

## Files

| File | Purpose |
|---|---|
| `MasterTBModel.R` | Main model — runs everything |
| `config.csv` | All 117 model parameters (values, distributions, sources) |
| `S2. Screening decision tree...xlsx` | Decision tree initial conditions (Zenner et al. 2025) |

## How to run

```bash
Rscript MasterTBModel.R
```

Or open in RStudio and run the full script. Runtime ~10–15 minutes (PSA: 1,000 simulations, parallelised).

Outputs are written to `output/` (PNG figures) and `output/csv/` (results tables).

## Model structure

- 16 health states, 43 screening strategies (34 sequential from decision tree + 8 parallel IGRA + no screening)
- 100,000 cohort, 660 monthly cycles (55-year lifetime horizon)
- NHS & PSS perspective; 3.5% discount rate (NICE reference case)
- PSA: 1,000 Monte Carlo simulations (beta/gamma distributions)
- DSA: tornado diagrams for 3 frontier strategies
- Background mortality: age-varying ONS life tables (25–79 years)

## Key outputs

| File | Contents |
|---|---|
| `output/csv/icer_table_basecase.csv` | Base-case ICERs for all 43 strategies |
| `output/csv/psa_results.csv` | PSA simulation results (cost, QALY per strategy per sim) |
| `output/csv/dsa_results.csv` | DSA results for tornado diagrams |
| `output/ceac.png` | Cost-effectiveness acceptability curve |
| `output/ce_plane.png` | Cost-effectiveness plane |

