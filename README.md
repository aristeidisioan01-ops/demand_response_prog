# Dynamic System Modeling & Optimization (MATLAB/Simulink)

This repository contains MATLAB scripts and Simulink models used for the modeling, simulation, and optimization of a dynamic energy system. Two different demand response (DR) schemes are implemented:

1. **Cost / CO₂ minimization:**  
   The Grey Wolf Optimizer (GWO) is used to minimize the objective function related to energy cost and CO₂ emissions.

2. **Loss minimization:**  
   The DR strategy aims at minimizing power losses in the network (power lines and transformer).

---

## Contents

- `codes/`  
  MATLAB scripts and functions.

- `models/`  
  Simulink models (`.slx`).

- `results/`  
  Simulation results (`.mat`, `.fig`).

- `docs/` (optional)  
  Additional documentation.

---

## Requirements

- MATLAB  
- Simulink  
- (Optional) Optimization / Global Optimization / Neural Network toolboxes, if required by the scripts.

---

## How to Run

### A. Base case – No demand response

1. Open MATLAB.
2. Run `initial_conditions.m`.
3. Open and run `no_demand_res.slx`.  
   After the simulation finishes, run `no_demand_results.m` to store the results.
4. Run `filter_no_signals.m` to extract the correct 24-hour values for each measurement.

### B. Demand response – Cost / CO₂ minimization

5. Run **one** of the following optimization scripts:
   - `gwo_nn_comb_co2.m` (CO₂-oriented objective), or  
   - `gwo_nn_comb_cost.m` (cost-oriented objective).
6. Open and run `with_demand_res_final_v5.slx`.  
   After the simulation finishes, run `with_demand_results.m`.
7. Run `filter_with_signals.m` to extract the correct 24-hour values for each measurement.  
   → Cost/CO₂-minimizing DR scenario completed.

### C. Demand response – Loss minimization (iterative GWO–simulation loop)

8. Close and reopen MATLAB (to clear the workspace).
9. Run `initial_conditions.m` to initialize the first iteration.
10. Open **only** `with_demand_res_final_v5_live.slx`.
11. Run `gwo_live_simulation_v2.m` and wait until all simulations and the GWO optimization are completed.
12. After completion, you will obtain:
    - `loss_history_figure.fig`
    - `load_compar.mat`
    - `load_history.mat`
13. Run `filter_with_loss_signals.m` to extract the correct 24-hour values for each measurement.
14. Run `calculate_cost_kgCO2.m` to compute the total cost and CO₂ emissions for this scenario.

---

## Comparison

15. Compare the two DR scenarios (cost/CO₂ minimization vs. loss minimization) using the generated `.mat` files and figures in `results/`.

