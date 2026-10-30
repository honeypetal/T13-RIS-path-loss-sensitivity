# How much does the direct-link path-loss model change simulated RIS gain?

MATLAB code and results for the paper *"How Much Does the Direct-Link Path-Loss Model Change Simulated RIS Gain? A Two-Model Sensitivity Study at 28 and 73 GHz"* by Aishwarya Vaidyanathan, Aanya Yadav and Shubham Kumar (School of Computer Science and Engineering, Vellore Institute of Technology).

## What is here

| File | Purpose |
|---|---|
| `scripts/run_revision.m` | **Main script for the paper.** Hook check, direct-link path-loss diagnostics (including the SimRIS street-canyon path-length test), main study with common random numbers (K = 1000), sensitivity sweep over the NLOS exponent, forced LOS/NLOS, equivalence tests, ablation and bandwidth check. Writes everything to `results_REV/`. |
| `scripts/geo_ris_channel.m` | Independent 3GPP TR 38.901-based geometric RIS channel model (same output format as SimRIS). |
| `scripts/poissrnd.m` | Poisson generator used by SimRIS when the Statistics Toolbox is not installed. |
| `scripts/ris_capacity.m`, `run_methodology1.m`, `run_methodology2.m`, `m2x_*.m`, `run_methodology3.m`, `run_extras.m` | Scripts from the earlier course-project stages (kept for reproducibility of the earlier report). |
| `results/results_REV/` | CSV files and figures produced by `run_revision.m`. |

## Requirements

* MATLAB R2024b (earlier releases from R2019b should work). No paid toolboxes are needed; `randsrc` from the Communications Toolbox is used by SimRIS itself.
* SimRIS Channel Simulator v2.0 by E. Basar and I. Yildirim: https://github.com/ifyildirim/SimRIS-Channel-Simulator (not redistributed here). Copy `SimRIS_GUI.m` and `SimRIS_GUI.fig` into the same folder as the scripts.

## How to reproduce

1. Put the files of `scripts/` and `SimRIS_GUI.m` / `SimRIS_GUI.fig` in one folder and make it MATLAB's current folder.
2. Open `run_revision.m`. Run once with `QUICK = true` (a few minutes) to check the set-up, then set `QUICK = false` and run again (about 20–60 minutes).
3. The CSV files and figures appear in `results_REV/`.

`run_revision.m` creates two files from `SimRIS_GUI.m` at run time: `SimRIS_ref.m` (an unmodified copy of SimRIS's channel function) and `SimRIS_r.m` (the same code with switches that are off unless set). With all switches off, both give identical channels for the same random seed; this is checked in part A of the script.

## Note on SimRIS's street-canyon direct link

In SimRIS v2.0 (`SimRIS_v18`, UMi branch of the direct-link generation), the length of each scattered direct-link path is computed as `a_c_rep_3 + norm(Tx_xyz - Coordinates2_3)`, i.e. the transmitter–scatterer distance twice, rather than transmitter–scatterer plus scatterer–receiver. `run_revision.m` reports results both for SimRIS as published and with this distance replaced by `norm(Rx_xyz - Coordinates2_3)` (switch `fix_umi`).

## Licence

The scripts in this repository are released under the MIT licence (see `LICENSE`). SimRIS is the work of its authors and is subject to its own terms.
