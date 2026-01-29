![geothermalwells_banner](https://github.com/user-attachments/assets/2f5a719c-245e-4476-88d6-8e972e5c5ff2)# GeothermalWells.jl

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cwittens.github.io/GeothermalWells.jl/stable/)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cwittens.github.io/GeothermalWells.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.12+-purple.svg)](https://julialang.org)
[![Build Status](https://github.com/cwittens/GeothermalWells.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cwittens/GeothermalWells.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cwittens/GeothermalWells.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cwittens/GeothermalWells.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

![Uplo<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 950 180" width="950" height="180">
  <defs>
    <style>
      @import url('https://fonts.googleapis.com/css2?family=Inter:wght@700&amp;display=swap');
      .title { font-family: 'Inter', 'Helvetica Neue', Arial, sans-serif; font-weight: 700; font-size: 72px; fill: #333; }
    </style>
  </defs>
  
  <!-- Logo scaled down, centered vertically -->
  <g transform="translate(10, 10) scale(0.49)">
    <!-- Connecting lines -->
    <line x1="175" y1="100" x2="88.398438" y2="250" stroke="#666" stroke-width="3" stroke-linecap="round"/>
    <line x1="175" y1="100" x2="261.601562" y2="250" stroke="#666" stroke-width="3" stroke-linecap="round"/>
    <line x1="88.398438" y1="250" x2="261.601562" y2="250" stroke="#666" stroke-width="3" stroke-linecap="round"/>
    
    <!-- Green well (TOP) -->
    <circle cx="175" cy="100" r="75" fill="#222"/>
    <circle cx="175" cy="100" r="73" fill="#389826"/>
    <circle cx="175" cy="100" r="62" fill="#2D7A5A"/>
    <circle cx="175" cy="100" r="38" fill="#389826"/>
    <circle cx="175" cy="100" r="25" fill="#6BC22F"/>
    
    <!-- Red well (BOTTOM LEFT) -->
    <circle cx="88.398438" cy="250" r="75" fill="#222"/>
    <circle cx="88.398438" cy="250" r="73" fill="#CB3C33"/>
    <circle cx="88.398438" cy="250" r="62" fill="#8B2525"/>
    <circle cx="88.398438" cy="250" r="38" fill="#CB3C33"/>
    <circle cx="88.398438" cy="250" r="25" fill="#F05A2A"/>
    
    <!-- Purple well (BOTTOM RIGHT) -->
    <circle cx="261.601562" cy="250" r="75" fill="#222"/>
    <circle cx="261.601562" cy="250" r="73" fill="#9558B2"/>
    <circle cx="261.601562" cy="250" r="62" fill="#5E4085"/>
    <circle cx="261.601562" cy="250" r="38" fill="#9558B2"/>
    <circle cx="261.601562" cy="250" r="25" fill="#B85FC5"/>
  </g>
  
  <!-- Text -->
  <text x="200" y="115" class="title">GeothermalWells.jl</text>
</svg>
ading geothermalwells_banner.svg…]()

**GeothermalWells.jl** is a Julia package for simulating deep borehole heat exchangers (DBHEs), supporting both single wells and well arrays. The package is designed to be easy to use, enabling rapid prototyping and well design exploration, while providing GPU-accelerated performance for long-term simulations.

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
  doi = {10.5281/zenodo.18405325}
}
```

## Authors

- [Collin Wittenstein](https://cwittens.github.io/) (Massachusetts Institute of Technology & Johannes Gutenberg University Mainz)

## License and Contributing

GeothermalWells.jl is licensed under the MIT license (see [LICENSE](LICENSE)).
Contributions are welcome - please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
