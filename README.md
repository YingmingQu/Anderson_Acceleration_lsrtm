# Anderson_Acceleration_lsrtm
# LSRTM with Anderson Acceleration (MATLAB)

## Description

This MATLAB program implements Least-Squares Reverse Time Migration (LSRTM) with Anderson Acceleration, a numerical optimization technique for seismic imaging. The algorithm improves convergence speed by leveraging historical iterations to accelerate updates, reducing computational costs in large-scale geophysical inverse problems. The code is designed for 2D seismic velocity model reconstruction using synthetic shot-gather data.

Key features:
- **Anderson Acceleration**: Enhances convergence of iterative gradient-based optimization.
- **Wavefield Modeling**: Supports acoustic wave propagation using finite-difference schemes.
- **Parallel Computing**: Utilizes MATLAB's `parfor` for shot-domain parallelization.

## Requirements
- MATLAB R2020a or later.
- Parallel Computing Toolbox (for shot-parallel processing).
- 8GB+ RAM (recommended for large models).

## Input Files
1. **`vpsmooth.dat`**: Smoothed velocity model (binary float32, size `[nz, nx]`).
2. **`ref.dat`**: True reflectivity model (binary float32, size `[nz, nx]`).
3. **`shot_born.dat`**: Synthetic shot-gather data (binary float32, size `[nt, nx * ns_hcp]`).

## Usage

### 1. Set Parameters
Modify the following parameters in `lsrtm_matlab_main.m`:
matlab
nx = 368;       % Grid size in x-direction

nz = 200;       % Grid size in z-direction

dx = 10.0;      % Grid spacing (meters)

dz = 10.0;      % Grid spacing (meters)

dt = 0.5;       % Time step (seconds)

tmax = 3.0;     % Maximum simulation time (seconds)

nit = 30;       % Number of iterations

frequency = 20; % Source peak frequency (Hz)

###2. Run the Code
Execute the script in MATLAB:

matlab
lsrtm_matlab_main; % Start LSRTM with Anderson Acceleration
###3. Outputs
mig_lsrtm[iter].dat: Reconstructed reflectivity model at each iteration.

gradient[iter].dat: Gradient updates.

residual.txt, misfit.txt: Convergence metrics.

Key Parameters
Parameter	Description
ns_hcp	Number of shots
crdist	Threshold distance for cutoffs
beta	Anderson acceleration mixing factor
mem_size	History size for acceleration

**Notes for Users**:
- Preprocess input velocity models using smoothing filters to avoid high-wavenumber artifacts.
- Adjust `mem_size` and `beta` to balance convergence speed and memory usage.
- Visualize outputs using MATLAB’s `imagesc` or external tools like Python’s `matplotlib`.

##Image
![Vp_model]()
![Shot]()
![Imaging](https://github.com/YingmingQu/Anderson_Acceleration_lsrtm/blob/main/mig_lsrtm30.jpg?raw=true)

