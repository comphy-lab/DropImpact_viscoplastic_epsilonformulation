# DropImpact_viscoplastic_epsilonformulation

This project contains C implementations for simulating two-phase flows, such as bouncing drops.

## Files in the Project

- `bounce.c`: Implementation of bouncing bubble simulation
- `two-phaseVP.h`: Header file for two-phase flow simulations
- `run_run.sh`: The script file for running the simulations in parallel. 
- `run_post.sh`: The script file for postprocess in parallel. 

## Prerequisites

To compile and run this project, you'll need:

- Basilisk C preinstalled. See: [Basilisk](https://basilisk.fr/)

# Required cases
## Phase Map Study
$We=10$, $We=4.5$ and $We=40$
```
tsnap=0.01
Ldomain=8.0
tmax=20.0
Bos=(0.0)
epsilons=(0.01)
MAXlevels=(9)
Wes=(10)
Ohs=(0.001 0.002 0.003 0.004 0.005 0.007 0.01 0.02 0.03 0.04 0.05 0.07 0.1 0.2 0.3 0.4 0.5 0.7 0.8 0.9 1.0 2.0) #22
Js=(0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2) #22
```


# Other tips
## viscoplastic droplet with bubble inside (May be interesting)
```
fraction(f, intersection(1. - R2Drop(x, y), R2Drop(x, y)-0.25));
```