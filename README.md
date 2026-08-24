# E-TANK replication package

This folder contains the standalone Dynare implementation used for the quantitative analysis in my thesis:

**Innovate or Compensate? Carbon-Tax Revenue Recycling in the Energy Transition**.

## Files

- `e_tank_model_public.mod` — standalone E-TANK perfect-foresight model.
- `README.md` — replication instructions.

## Requirements

- MATLAB
- Dynare (the public file was validated with Dynare 6.5)

## Horizons

- Main economic reporting horizon: **30 years = 120 quarters**.
- Perfect-foresight numerical solution horizon: **400 years = 1,600 quarters**.

The 400-year horizon is used to minimize terminal-horizon effects and should not be interpreted as a 400-year economic forecast.

## Reference outputs

The run reproduces the final thesis results:

| Outcome | Reference value |
|---|---:|
| Year-30 domestic emissions reduction | 36.940975% |
| Year-30 net-output deviation | -8.614750% |
| Year-30 clean knowledge, `A_c` | 10.482160 |
| Year-30 clean-energy share | 40.497866% |
| HtM lifetime S-CEV | 114.118122% |
| Ricardian lifetime S-CEV | -18.629561% |
| Aggregate lifetime money metric | 5.801250 |

The initial HtM energy burden is 0.18 and the reference per-capita transfer ratio `T_H/T_R` is 16.

## Model scope

E-TANK is a two-agent environmental real general-equilibrium model with hand-to-mouth and Ricardian households, Stone-Geary energy demand, clean and dirty energy, endogenous clean-knowledge accumulation, nested-CES production, carbon-tax revenue recycling, a small-open-economy trade closure, and reduced-form climate damages.
