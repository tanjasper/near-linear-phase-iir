# near-linear-phase-iir

This repo contains code for the following paper:  
J. Tan and C. S. Burrus, “Near-linear-phase IIR filters using Gauss-Newton optimization,” in IEEE Int. Midwest Symp. Circuits Syst., Aug. 2019  
Paper link: https://ieeexplore.ieee.org/abstract/document/8885116

The function `gauss_newton_iir.m` is our implementation of Algorithm 1 in the paper. It generates an IIR filter using the Gauss-Newton method.  
The `generate_paper_results.m` script reproduces the results from the paper (specifically, Fig. 1, Table I, and Fig. 2).
This script requires the generation of multiple IIR filters of different configurations, which can take hours.
Thus, these filters have been saved beforehand in `data/filters_diff_zeros_poles_L1024_fp0p125_fs0p135_iters3000.mat`.
The script to generate these filters is `generate_paper_filters.m`.
Note that running `generate_paper_filters.m` may take a few hours and is currently set to overwrite the pre-saved filters in the `data` folder.

This code has been tested on Matlab R2018a running on Ubuntu 16.04.

Contact: Jasper Tan (jaspertan@rice.edu)
