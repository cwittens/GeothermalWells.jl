# GeothermalWells.jl

```@raw html
<p align="center">
  <img height="200px" alt="logo" src="assets/geothermalwells_banner.svg">
</p>
```

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cwittens.github.io/GeothermalWells.jl/stable/)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cwittens.github.io/GeothermalWells.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.12+-purple.svg)](https://julialang.org)
[![Build Status](https://github.com/cwittens/GeothermalWells.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cwittens/GeothermalWells.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cwittens/GeothermalWells.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cwittens/GeothermalWells.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

**GeothermalWells.jl** is a Julia package for simulating deep borehole heat exchangers (DBHEs), supporting both single wells and well arrays. The package is designed to be easy to use, enabling rapid prototyping and well design exploration, while providing GPU-accelerated performance for long-term simulations.

!!! note "Under Development"
    This package is currently under active development and there may be breaking changes.

## Features

- Full 3D simulation of coaxial deep borehole heat exchangers
- Support for single wells and multi-well arrays
- GPU acceleration via [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) (CUDA, AMD ROCm)
- Adaptive grid generation with fine resolution near boreholes
- Integration with [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) solvers

## Installation

GeothermalWells.jl works with Julia v1.12 and newer. Install it from the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add("GeothermalWells")
```

For GPU support, also install the appropriate backend:
```julia
julia> Pkg.add("CUDA")        # for NVIDIA GPUs
julia> Pkg.add("AMDGPU")      # for AMD GPUs
```

## Quick Example

```julia
using GeothermalWells
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
using KernelAbstractions: CPU

# Define material properties
materials = HomogenousMaterialProperties{Float64}(
    2.88,    # k_rock [W/(m·K)]
    2.17e6,  # rho_c_rock [J/(m³·K)]
    0.6,     # k_water [W/(m·K)]
    4.18e6,  # rho_c_water [J/(m³·K)]
    44.5,    # k_steel [W/(m·K)]
    3.73e6,  # rho_c_steel [J/(m³·K)]
    0.26,    # k_insulating [W/(m·K)]
    1.96e6,  # rho_c_insulating [J/(m³·K)]
    1.0,     # k_backfill [W/(m·K)]
    1.0      # rho_c_backfill [J/(m³·K)]
)

# Define borehole geometry
borehole = Borehole{Float64}(
    0.0, 0.0,    # center position (x, y) [m]
    2000.0,      # depth [m]
    0.0381,      # inner pipe radius [m]
    0.01,        # inner pipe thickness [m]
    0.0889,      # outer pipe radius [m]
    0.01,        # outer pipe thickness [m]
    0.0989,      # backfill radius [m]
    10.0,        # mass flow rate [kg/s]
    500.0        # insulation depth [m]
)

# Create grids
gridx = create_adaptive_grid_1d(xmin=-100.0, xmax=100.0, dx_fine=0.0025,
    growth_factor=1.3, dx_max=10.0, boreholes=(borehole,),
    backend=CPU(), Float_used=Float64, direction=:x)
gridy = create_adaptive_grid_1d(xmin=-100.0, xmax=100.0, dx_fine=0.0025,
    growth_factor=1.3, dx_max=10.0, boreholes=(borehole,),
    backend=CPU(), Float_used=Float64, direction=:y)
gridz = create_uniform_gridz_with_borehole_depths(zmin=0.0, zmax=2200.0,
    dz=100.0, boreholes=(borehole,), backend=CPU())

# Initial temperature with geothermal gradient
T0 = initial_condition_thermal_gradient(CPU(), Float64, gridx, gridy, gridz;
    T_surface=10.0, gradient=0.03)

# Create cache and solve
cache = create_cache(backend=CPU(), gridx=gridx, gridy=gridy, gridz=gridz,
    materials=materials, boreholes=(borehole,),
    inlet_model=ConstantInlet{Float64}(20.0))

callback, saved_values = get_simulation_callback(saveat=[0.0, 3600.0])

prob = ODEProblem(rhs_diffusion_z!, T0, (0.0, 3600.0), cache)
solve(prob, ROCK2(eigen_est=eigen_estimator), dt=60.0, callback=callback, adaptive=false)
```

See the [Getting Started](@ref getting-started) guide for a complete walkthrough.

## Examples and Reproducibility

For more comprehensive examples, including validation against published literature, see the
[reproducibility repository](https://github.com/cwittens/2026_DBHEs_Arrays). This repository contains:

- Reproduction of results from Li et al. (2021), Hu et al. (2020), Brown et al. (2023)
- Convergence studies
- Well array simulations

## Referencing

If you use [GeothermalWells.jl](https://github.com/cwittens/GeothermalWells.jl) in your research, please cite it:

```bibtex
@software{wittenstein2026geothermalwells,
  title={{GeothermalWells.jl}: {GPU}-accelerated simulation of deep borehole heat exchanger (DBHE) arrays},
  author={Wittenstein, Collin},
  year={2026},
  howpublished={\url{https://github.com/cwittens/GeothermalWells.jl}},
  doi = {10.5281/zenodo.18405325}
}
```

## Authors

- [Collin Wittenstein](https://cwittens.github.io/) (Massachusetts Institute of Technology & Johannes Gutenberg University Mainz)

## License

GeothermalWells.jl is licensed under the MIT license (see [License](@ref license-page)).
