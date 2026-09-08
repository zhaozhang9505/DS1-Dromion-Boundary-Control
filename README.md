# DS1-Dromion-Boundary-Control

This repository contains MATLAB programs for boundary-driven circular transport of DS-I dromions.

- `circular_Dromion.m` simulates one period of circular motion from an exact initial condition and compares the numerical intensity with the analytical result.
- `gaussian_Dromion_smooth_damped.m` simulates ten periods from a Gaussian initial condition with smooth boundary damping.
- `plot_gaussian_Dromion_diagnostics.m` generates plots from the saved Gaussian simulation data without repeating the simulation.

The programs use Fourier pseudospectral differentiation and fourth-order Runge–Kutta time stepping.
Place all three files in the same folder and set it as the MATLAB current folder.
Run either main script; the Gaussian script calls the plotting function automatically.
Results are saved in automatically created subfolders, and repeated runs may overwrite existing output files.
