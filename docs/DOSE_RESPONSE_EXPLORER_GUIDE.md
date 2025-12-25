# Dose-Response Explorer Guide

## Overview
The **Dose-Response Explorer** helps you fit potency curves and analyze longitudinal readouts directly in the dashboard. It supports classical **3-parameter logistic (3PL)** and **4-parameter logistic (4PL)** fits (SciPy), plus an optional **Bayesian 4PL** fit (PyMC) for cases where you want credible intervals and more robust uncertainty handling.

You can use it in four modes: **Single Curve**, **Compare Curves**, **Time Series**, and **Batch Analysis**.

## Single Curve Fitting
### Data input
Provide **(concentration, response)** pairs by:
- **Manual entry**: one pair per line as `concentration,response` (commas, spaces, and tabs are accepted).
- **CSV upload**: a file with columns `concentration` and `response`.

### Model selection
Choose one of:
- **3PL**: estimates `top`, `ec50`, `hill_slope` (no `bottom` parameter).
- **4PL**: estimates `bottom`, `top`, `ec50`, `hill_slope`.
- **Bayesian 4PL**: optional dependency; returns `ec50` mean + credible interval and `hill_slope` (MVP response does not include full top/bottom summaries).

### Interpreting results
- **EC50**: concentration at half-maximal effect (lower EC50 = higher potency).
- **Hill slope**: steepness of the transition.
- **Top / Bottom**: response asymptotes (4PL has both; 3PL assumes bottom is ~0 in this simplified form).
- **R²**: goodness-of-fit for frequentist fits (not returned for Bayesian in this MVP).

The chart uses a **log x-axis** and plots raw points plus the fitted curve. If confidence bands are available, they appear as a shaded region.

## Comparing Curves
To compare multiple compounds/conditions:
- Upload multiple CSV files (each becomes one curve), or
- Provide manual triplets `compound_id, concentration, response`.

The comparison view overlays fitted curves with distinct colors and provides a table for **EC50 / Hill / R²** across curves to spot potency shifts.

## Time Series Analysis
The Time Series tab takes **(timepoint, value)** pairs from manual input or CSV (`timepoint`, `value`) and reports:
- **Slope** and **p-value** from a linear trend test
- **Direction**: increasing / decreasing / stable

Optional features:
- **Smoothing**:
  - Savitzky–Golay (SciPy) for light denoising
  - LOWESS (statsmodels, optional) for nonparametric smoothing
- **Changepoints**: a simple rolling-variance heuristic that highlights time indices where variability spikes.

## Batch Processing
Batch fitting expects a CSV with:
- `compound_id`
- `concentration_nm`
- `response_percent`

The dashboard groups rows by `compound_id` and fits up to 50 curves per run. It then displays:
- A summary table (including warnings for problematic fits)
- A parameter heatmap (EC50 / Hill / R²) for quick triage

## API Reference
All endpoints are under `/api/explorer`:

### 1) Fit a single dose-response curve
`POST /api/explorer/dose-response/fit`

Example request:
```json
{
  "concentrations": [1, 10, 100, 1000],
  "responses": [0.1, 0.2, 0.8, 0.95],
  "model": "4PL"
}
```

### 2) Fit and compare multiple curves
`POST /api/explorer/dose-response/compare`

Example request:
```json
{
  "fits": [
    {"concentrations": [1, 10, 100, 1000], "responses": [0.1, 0.2, 0.8, 0.95], "model": "4PL"},
    {"concentrations": [1, 10, 100, 1000], "responses": [0.2, 0.3, 0.7, 0.9], "model": "4PL"}
  ]
}
```

### 3) Analyze a time series
`POST /api/explorer/timeseries/analyze`

Example request:
```json
{
  "values": [1, 2, 3],
  "timepoints": [0, 1, 2],
  "smooth_method": "savgol",
  "detect_changepoints": true,
  "changepoint_threshold": 2.0
}
```

### 4) Compare trajectories
`POST /api/explorer/timeseries/compare`

Example request:
```json
{
  "series": [
    {"values": [1, 2, 3], "timepoints": [0, 1, 2]},
    {"values": [3, 2, 1], "timepoints": [0, 1, 2]}
  ]
}
```

## Tips
- Use **≥ 6 dose points** when possible; sparse curves can fit poorly and may produce warnings.
- Make sure **concentrations are > 0** (logistic models assume positive concentrations).
- Treat EC50 estimates with caution if the curve does not visibly reach a plateau (top/bottom poorly constrained).


