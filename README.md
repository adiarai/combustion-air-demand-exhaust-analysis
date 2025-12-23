 # Energy Process and Systems Engineering – Assignment 2

## 📘 Project Description
This repository contains Assignment 2 for the course *Energy Process and Systems Engineering*.  
The project focuses on stoichiometric combustion analysis, air demand calculation, and exhaust gas composition using both molecular and elemental approaches.

## 🔹 Problem Overview

### Problem 1: Air Quantity – Mole Fractions (Molecular Mode)
- Calculation of minimum and actual air demand
- Based on mole fractions of fuel components
- Includes excess air factor (λ)

### Problem 2: Air Quantity – Mass Fractions (Elemental Mode)
- Air demand calculation using elemental mass fractions (C, H, S, O)
- Conversion between mass and mole fractions
- Element-based combustion reactions

### Problem 3: Exhaust Gas – Mole Fractions
- Molar balance using stoichiometric matrices
- Exhaust gas composition based on extent of reaction

### Problem 4: Exhaust Gas – Mass Fractions
- Mass balance approach for exhaust gas
- Conversion between mass and mole fractions

## 🧮 Implementation
The project is implemented in **Julia** and consists of:
- `fun.jl`: Core functions for stoichiometry, air demand, and exhaust gas calculations
- `main.jl`: Execution scripts for natural gas and coal combustion, including plots

## 🛠 Tools & Methods
- Julia programming language
- Linear algebra and stoichiometric matrices
- Mole and mass balance calculations
- Visualization using Plots.jl

## 📈 Output
- Air demand vs. excess air ratio (λ)
- Exhaust gas composition (molar and mass fractions)

## 📌 Academic Context
This repository is intended for academic use and documentation of coursework.
## 🧮 Actual Numerical Results

### Natural Gas Combustion (Molar Basis)

At stoichiometric conditions (λ = 1.0):

- Minimum air demand:  
  X_air,min = 9.52 mol air / mol fuel

- Oxygen supplied:  
  X_O₂,in = 2.00 mol O₂ / mol fuel

- Exhaust gas mole fractions:
  - O₂ ≈ 0.00
  - N₂ ≈ 0.71
  - CO₂ ≈ 0.09
  - H₂O ≈ 0.20

At excess air ratio λ = 1.5:

- Actual air demand:  
  X_air,in = 14.28 mol air / mol fuel

- Exhaust gas mole fractions:
  - O₂ ≈ 0.07
  - N₂ ≈ 0.74
  - CO₂ ≈ 0.06
  - H₂O ≈ 0.13
### Coal Combustion (Mass Basis)

At λ = 1.0:

- Minimum air demand:  
  W_air,min = 8.64 kg air / kg fuel

- Exhaust gas mass fractions:
  - CO₂ ≈ 0.18
  - H₂O ≈ 0.11
  - SO₂ ≈ 0.02
  - N₂ ≈ 0.69

At λ = 1.8:

- Actual air demand:  
  W_air,in = 15.55 kg air / kg fuel

- Residual O₂ mass fraction ≈ 0.08

