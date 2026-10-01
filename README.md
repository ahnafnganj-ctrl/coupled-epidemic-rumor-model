# Coupled Epidemic–Rumor Model: MATLAB Codes

This repository contains MATLAB codes for numerical simulation, threshold analysis, case studies, heatmap visualization, and sensitivity analysis of a deterministic coupled epidemic–rumor model with imperfect vaccination and behavioral feedback.

## Repository Structure

```
Coupled-Epidemic-Rumor-Model/
│
├── README.md
│
├── baseline_trajectories.m
├── Case_studies.m
├── coupling_comparison.m
├── Long_term_infection_rumor_transmission.m
├── Heatmap_analysis.m
├── PRCC_analysis.m
│
└── matlab_output/
    └── data/
        ├── trajectory_*.csv
        ├── response_*.csv
        └── map_*.csv
```

## MATLAB Files

### `baseline_trajectories.m`
Generates baseline epidemic and rumor trajectories, including extinction, persistence, and coexistence scenarios.

### `Case_studies.m`
Simulates different behavioral and epidemiological scenarios:
- Baseline
- Low rumor transmission
- High rumor transmission
- High vaccination uptake
- High vaccine leakiness

### `coupling_comparison.m`
Compares the effects of epidemic–rumor coupling mechanisms and evaluates their impact on infection dynamics.

### `Long_term_infection_rumor_transmission.m`
Generates long-term parameter-response analyses and threshold-based transmission figures.

### `Heatmap_analysis.m`
Produces two-parameter heatmaps showing the effects of parameter interactions on final infection, vaccination, and rumor levels.

### `PRCC_analysis.m`
Performs global sensitivity analysis using Latin Hypercube Sampling (LHS) and Partial Rank Correlation Coefficient (PRCC) for:

- \(R_{0E}\): epidemic reproduction number
- \(R_{0R}\): rumor reproduction number
- \(R_{E|1}\): disease invasion number around the persistent-rumor equilibrium

## Requirements

- MATLAB R2020a or later
- Standard MATLAB functions only

## Running the Codes

Run individual scripts directly in MATLAB:

```matlab
baseline_trajectories
```

```matlab
Case_studies
```

```matlab
coupling_comparison
```

```matlab
Long_term_infection_rumor_transmission
```

```matlab
Heatmap_analysis
```

```matlab
results = PRCC_analysis;
```

## Data

The visualization scripts require CSV files located in:

```
data
```

including trajectory, response, and heatmap datasets.

## Citation

If you use these codes, please cite the associated manuscript:

**"Dynamics of Epidemics and Rumors: A Coupled System with Imperfect Vaccination and Behavioral Feedback."**
