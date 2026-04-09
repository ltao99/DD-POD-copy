# Accelerating 3D Magnetotelluric Forward Modelling with Domain Decomposition and Order-Reduction Methods

**Authors/Information:** Luis Tao and Fabio Zyserman  
**Paper reference:** ** (https://doi.org/10.31223/X5JV1C)*

## Abstract
This repository contains the source code, sample setups, and data for the 3D Magnetotelluric (MT) forward modelling code utilizing Domain Decomposition (DD) and Proper Orthogonal Decomposition (POD) techniques. The purpose of this software is to speed up complex 3D MT simulations by reducing the order of the problem while preserving acceptable accuracy in apparent resistivity and phase. 

## Technical Context: Domain Decomposition & POD
The computational domain of the 3D MT forward modelling is computationally expensive, forming a large-scale linear system. To efficiently solve it the following methods are employed:

1. **Domain Decomposition (DD):** The code divides the global 3D computational domain into several local subdomains which are solved independently and coupled iteratively at the interfaces. Specifically, this domain decomposition strategy partitions the domain **only in the horizontal directions**, meaning that the subdomains only share horizontal interfaces. This distributes the memory footprint and the computational load across multiple processors (using MPI). One processor is assigned per subdomain and polarization.
2. **Proper Orthogonal Decomposition (POD / Order-Reduction):** In addition to the full-order formulation, the DD solver can utilize Reduced Basis-Proper Orthogonal Decomposition (RB-POD) to accelerate computation. By extracting POD bases from pre-computed high-fidelity snapshots and projecting the solver components into a lower-dimensional space, the scheme drastically decreases the Degrees of Freedom (DoF).

## Inputs

The DD solver requires input files and parameters defining the physical model and numerical setup.

### 1. `sizes.f90` (Model & Solver Configuration)
This file controls the domain partitioning and basic mesh dimensions.
* `nsy, nsz`: Number of subdomains in each direction. (`nsy` is 1).
* `nprocs_input`: Total processors used. Must be `2 × nsy × nsz` (two processors are required per subdomain, one for each polarization).
* `ngx, ngy, ngz`: Number of elements along the X, Y, Z directions. `ngz` must be a multiple of `nsz`.

### 2. `input_mesh_MT` (Discretization & Frequencies)
Defines the spatial discretization and frequency sampling. The file structure starts with the number of frequencies and their values. Afterward, for each dimension (X, Y, Z), it defines coordinates of the macro-nodes in meters and the number of partitions between them.

### 3. `conduc_input` (Conductivity Model)
Defines the conductivity distribution (in S/m) of each element in the 3D mesh (X × Y × Z dimensions). The file appends a 1D vertically-layered background model, consisting of layer thicknesses (where the last layer has no thickness assigned) and layer conductivities. Air layers must be explicitly included.

### 4. `MOD_Constants.f90` (Physical Constants & Hyperparameters)
Controls full-order or reduced-order flags and DD hyper-parameters:
* `usepod`: Set to `.TRUE.` for Reduced-Order Modelling, `.FALSE.` for Full-Order.
* **POD settings**: `Nsnapshots` (number of snapshots), `tolsvd` (SVD truncation tolerance), `podpath` (directory for snapshot matrices folder.
* **DD settings**: `relx` (relaxation parameter), `penal` (penalization factor), `maxiter`, and `convergence_tol`.

## Outputs

The solver evaluates the physical quantities and dumps them to the outputs directory:
1. **Electromagnetic Fields:** Electric and magnetic field components evaluated at the midpoints of the surface elements stored in the /outputs.
2. **Apparent Resistivity & Phase:** Computed for each polarization evaluated at the midpoints of the surface mesh stored in /outputs.
3. **Performance Logs:** Runtime metrics, timings for each computational phase, and iterative convergence progress stored in /outputs.

## Repository Structure & Data
* **`src/`**: Fortran source code `*.f90`.
* **`examples/`**: Contains configurations (`conduc_input`, `input_mesh_MT`), scripts (`job_mt_monolithic.sh`), and standard parameter setups for 20x20x18, 40x40x36, and 80x80x72 mesh executions. Those example folders define the inputs to generate full-order models.
* **`data/`**: Stores the snapshots databases. `data/20x20x18/Snapshots` contains the high-fidelity state snapshots required to deploy the POD-based reduction.
* **`docs/`**: Includes the full PDF `Code_Guide.pdf` for extended context.

## Prerequisites
- **Fortran Compiler:** Intel `mpif90` / Intel oneAPI
- **MPI:** For distributed Domain Decomposition
- **MKL (Math Kernel Library):** For BLAS and LAPACK routines
- **MUMPS & PARMETIS/METIS:** For the distributed sparse direct solver

## Installation & Usage
1. Configure build environment: Review `Makefile` paths (`PATHLIB`, `MKLROOT`, `PATHCOMPILER`) to match your compiler deployment.
2. Build the program: `make`
3. Execute via MPI. For instance, `nsy=1` and `nsz=3` uses 6 processes:
```bash
mpirun -np 6 ./max3d-par
```
*(Ensure that all required input files the code targets, such as `conduc_input`, `input_mesh_MT`, and others are copied to the working directory alongside the executable).*

## Acknowledgments
Development was carried out as part of the Marie Skłodowska-Curie EarthSafe project (EU Horizon Europe programme grant No. 101120556). Special thanks to Prof. Fabio Zyserman for the foundation code and Dr. Constanza Manassero for stiffness matrix assemblage routines.
