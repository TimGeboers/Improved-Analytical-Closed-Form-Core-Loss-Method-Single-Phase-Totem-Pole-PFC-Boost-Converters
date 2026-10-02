# Improved Analytical Closed-Form Core Loss Calculation Method for Single-Phase Totem-Pole PFC Boost Converters

This repository contains the MATLAB implementation for calculating core losses in coupled inductors for single-phase totem-pole PFC boost converters using analytical closed-form improved Generalized Steinmetz Equations (iGSE).

## Citation
T. Geboers, W. Vanderwegen, W. Martinez and C. Suarez, "Improved Core Loss Calculation Method for Single-Phase Totem-Pole PFC Boost Converters," in IEEE Energy Conversion Congress and Exposition (ECCE), 2026, doi: **not yet available**

## Contributors
- **Tim Geboers**  
  [tim.geboers1@kuleuven.be](mailto:tim.geboers1@kuleuven.be) — KU Leuven / EnergyVille
- **Wout Vanderwegen**  
  [wout.vanderwegen@kuleuven.be](mailto:wout.vanderwegen@kuleuven.be) — KU Leuven / EnergyVille
- **Wilmar Martinez**  
  [wilmar.martinez@kuleuven.be](mailto:wilmar.martinez@kuleuven.be) — KU Leuven / EnergyVille
- **Camilo Suarez**  
  [camilo.suarez@kuleuven.be](mailto:camilo.suarez@kuleuven.be) — KU Leuven / EnergyVille

## Contents
- Analytical closed-form iGSE expressions separating major and minor B-H loop core losses
- Analytical integral iGSE baseline calculation routine for numerical verification

## 🔓 Read the Paper
- On [ResearchGate](https://www.researchgate.net/publication/414649223) (Open Access pre-print)
- On [IEEE Xplore] () (Coming Soon)

## Repository Structure
- `/code/*.m files`: Analysis and core calculation script
  - [`Improved_Analytical_Closed_Form_Core_Loss_Method_Single_Phase_Totem_Pole_PFC_Boost_Converters.m`](Improved_Analytical_Closed_Form_Core_Loss_Method_Single_Phase_Totem_Pole_PFC_Boost_Converters.m) Main calculation script for closed-form and integral iGSE core loss models

## Requirements
- MATLAB R2025b

### Individual Script Execution
1. Clone repository
2. Run any `.m` file in MATLAB R2025b
   - Files follow paper's sequence
   - Executable independently
