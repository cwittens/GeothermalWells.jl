# GeothermalWells.jl

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cwittens.github.io/GeothermalWells.jl/stable/)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cwittens.github.io/GeothermalWells.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.12+-purple.svg)](https://julialang.org)
[![Build Status](https://github.com/cwittens/GeothermalWells.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cwittens/GeothermalWells.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cwittens/GeothermalWells.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cwittens/GeothermalWells.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

**GeothermalWells.jl** is a Julia package for simulating Deep Borehole Heat Exchanger (DBHE) Arrays.

> **Note:** This package is currently under active development and there may be breaking changes.

## Installation

GeothermalWells.jl works with Julia v1.12 and newer. Install it from the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add("GeothermalWells")
```

## Quick Start (pseudocode)

```julia
using GeothermalWells
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
using KernelAbstractions: CPU

# Set up material properties and borehole geometry
materials = HomogenousMaterialProperties{Float64}(...)
borehole = Borehole{Float64}(...)

# Create adaptive grids
gridx = create_adaptive_grid_1d(...)
gridy = create_adaptive_grid_1d(...)
gridz = create_uniform_gridz_with_borehole_depths(...)

# Initial temperature field with geothermal gradient
T0 = initial_condition_thermal_gradient(backend, Float64, gridx, gridy, gridz;
    T_surface=2.29, gradient=0.035)

# Create simulation and solve
cache = create_cache(backend=CPU(), gridx=gridx, gridy=gridy, gridz=gridz,
    materials=materials, boreholes=(borehole,), inlet_model=ConstantInlet{Float64}(20.0))
prob = ODEProblem(rhs_diffusion_z!, T0, (0.0, 3600.0), cache)
callback, saved_values = get_simulation_callback(...)
solve(prob, ROCK2(), dt=60.0, callback=callback)
```

See the [examples/](examples/) folder for complete working examples.

## Documentation

For more details, see the [documentation](https://cwittens.github.io/GeothermalWells.jl/stable/).

If you use [GeothermalWells.jl](https://github.com/cwittens/GeothermalWells.jl) in your research, please cite it:

```bibtex
@software{wittenstein2026geothermalwells,
  title={{GeothermalWells.jl}: {GPU}-accelerated simulation of deep borehole heat exchanger (DBHE) arrays},
  author={Wittenstein, Collin},
  year={2026},
  howpublished={\url{https://github.com/cwittens/GeothermalWells.jl}},
  doi = {10.5281/zenodo.YYYYYYY}
}
```

TODO Zenodo.

## Authors

- [Collin Wittenstein](https://cwittens.github.io/) (Massachusetts Institute of Technology & Johannes Gutenberg University Mainz)

## License and Contributing

GeothermalWells.jl is licensed under the MIT license (see [LICENSE](LICENSE)).
Contributions are welcome - please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
