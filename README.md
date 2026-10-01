# Dynamics of Epidemics and Rumors: MATLAB Simulation and Sensitivity Analysis Codes

This repository contains MATLAB codes for numerical simulation, equilibrium analysis, threshold analysis, parameter sensitivity analysis, and visualization of a deterministic coupled epidemic–rumor model with imperfect vaccination and behavioral feedback.

The model couples disease transmission dynamics with information/rumor propagation. The epidemic subsystem includes susceptible, infected, vaccinated, and recovered populations, while the rumor subsystem describes unaware, active rumor spreaders, and corrected individuals.

## Model Framework

The coupled system considers:

- Disease states:
  - S: susceptible population
  - I: infected population
  - V: vaccinated population
  - R: recovered population

- Rumor states:
  - U: unaware population
  - A: active rumor spreaders
  - C: corrected/informed population

Key model features:

- Imperfect vaccination
- Rumor-induced reduction of vaccination uptake
- Infection-driven rumor activation
- Bidirectional epidemic–information feedback

All simulations use normalized model time and illustrative parameter values described in the associated manuscript.

---

# MATLAB Requirements

Recommended:

- MATLAB R2020a or later
- Base MATLAB functions only

Required functionality:

- ODE45 solver
- Standard plotting functions
- Linear algebra functions

No additional MATLAB toolbox is required for the main simulations.

---

# Repository Structure
