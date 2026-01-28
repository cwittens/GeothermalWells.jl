# [Getting Started](@id getting-started)

This guide walks through a complete simulation of a single deep borehole heat exchanger.

## Overview

A typical GeothermalWells.jl simulation involves:
1. Defining material properties (rock, water, steel, etc.)
2. Specifying borehole geometry
3. Creating computational grids
4. Setting initial conditions
5. Choosing an inlet model
6. Running the simulation

## Complete Example

```julia
using GeothermalWells
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
using KernelAbstractions: CPU

# Choose backend: CPU() for testing, CUDABackend() for NVIDIA GPU
backend = CPU()
Float_used = Float64
```

!!! warning "Use Float64"
    Using `Float32` can cause numerical instabilities in some cases. It is recommended to use `Float64` for reliable results.

### Material Properties

Define thermal properties for all materials in the simulation:

```julia
materials = HomogenousMaterialProperties{Float_used}(
    2.88,       # k_rock - rock thermal conductivity [W/(m·K)]
    2.17e6,     # rho_c_rock - rock volumetric heat capacity [J/(m³·K)]
    0.6,        # k_water [W/(m·K)]
    4.18e6,     # rho_c_water [J/(m³·K)]
    44.5,       # k_steel [W/(m·K)]
    3.73e6,     # rho_c_steel [J/(m³·K)]
    0.26,       # k_insulating [W/(m·K)]
    1.96e6,     # rho_c_insulating [J/(m³·K)]
    1.0,        # k_backfill [W/(m·K)]
    1.0         # rho_c_backfill [J/(m³·K)]
)
```

For layered rock formations, use [`StratifiedMaterialProperties`](@ref) instead.

### Borehole Geometry

Define the coaxial borehole heat exchanger geometry:

```julia
borehole = Borehole{Float_used}(
    0.0,             # xc - x-coordinate of center [m]
    0.0,             # yc - y-coordinate of center [m]
    2000.0,          # h - borehole depth [m]
    0.0381,          # r_inner - inner pipe radius [m]
    0.01,            # t_inner - inner pipe wall thickness [m]
    0.0889,          # r_outer - outer pipe inner radius [m]
    0.01,            # t_outer - outer pipe wall thickness [m]
    0.0989,          # r_backfill - backfill outer radius [m]
    10.0,            # ṁ - mass flow rate [kg/s]
    500.0            # insulation_depth - depth of inner pipe insulation [m]
)

boreholes = (borehole,)  # tuple of all boreholes
```

### Grid Setup

Create adaptive grids that are fine near the borehole and coarse far away:

```julia
# Domain boundaries
xmin, xmax = -100.0, 100.0
ymin, ymax = -100.0, 100.0
zmin, zmax = 0.0, 2200.0

# Grid parameters
dx_fine = 0.0025      # fine spacing near borehole [m]
growth_factor = 1.3   # geometric growth rate
dx_max = 10.0         # maximum spacing far from borehole [m]
dz = 100.0            # vertical spacing [m]

gridx = create_adaptive_grid_1d(
    xmin=xmin, xmax=xmax,
    dx_fine=dx_fine, growth_factor=growth_factor, dx_max=dx_max,
    boreholes=boreholes, backend=backend, Float_used=Float_used, direction=:x
)

gridy = create_adaptive_grid_1d(
    xmin=ymin, xmax=ymax,
    dx_fine=dx_fine, growth_factor=growth_factor, dx_max=dx_max,
    boreholes=boreholes, backend=backend, Float_used=Float_used, direction=:y
)

gridz = create_uniform_gridz_with_borehole_depths(
    zmin=zmin, zmax=zmax, dz=dz,
    boreholes=boreholes, backend=backend
)

println("Grid size: $(length(gridx)) x $(length(gridy)) x $(length(gridz))")
```

### Initial Condition

Set the initial temperature field with a geothermal gradient:

```julia
T0 = initial_condition_thermal_gradient(
    backend, Float_used, gridx, gridy, gridz;
    T_surface=10.0,    # surface temperature [°C]
    gradient=0.03      # thermal gradient [K/m]
)
```

### Inlet Model

Choose how the inlet temperature is determined:

```julia
# Option 1: Constant inlet temperature
inlet_model = ConstantInlet{Float_used}(20.0)  # 20°C

# Option 2: Heat exchanger (inlet = outlet - ΔT)
# inlet_model = HeatExchangerInlet{Float_used}(5.0)  # ΔT = 5K
```

### Create Cache and Solve

```julia
cache = create_cache(
    backend=backend,
    gridx=gridx,
    gridy=gridy,
    gridz=gridz,
    materials=materials,
    boreholes=boreholes,
    inlet_model=inlet_model
)

# Time span: simulate for 1 day
tspan = (0.0, 3600.0 * 24.0)

prob = ODEProblem(rhs_diffusion_z!, T0, tspan, cache)

# Callback handles ADI (horizontal diffusion) and advection, plus saves solutions
# This is required - without it, only vertical diffusion is computed
callback, saved_values = get_simulation_callback(
    saveat=range(tspan..., 10),  # times to save solutions
    print_every_n=100            # print progress every N steps
)

# Solve with ROCK2 (stabilized explicit Runge-Kutta)
Δt = 80.0  # time step [s]

solve(
    prob,
    ROCK2(max_stages=100, eigen_est=eigen_estimator),
    save_everystep=false,
    callback=callback,
    adaptive=false,
    dt=Δt,
    maxiters=Int(1e10)
)
```

## GPU Acceleration

For GPU acceleration, simply change the backend:

```julia
using CUDA: CUDABackend

backend = CUDABackend()

# All other code remains the same!
```

## Well Arrays

To simulate multiple boreholes, create a tuple of [`Borehole`](@ref) objects:

```julia
# 3x3 array with 30m spacing
boreholes = tuple(
    (Borehole{Float_used}(
        xc, yc,          # center coordinates
        2000.0,          # depth
        # ... other parameters
    ) for xc in [-30.0, 0.0, 30.0], yc in [-30.0, 0.0, 30.0])...
)
```

## Next Steps

- See the [API Reference](@ref api-reference) for detailed documentation
- Check the `examples/` folder:
  - `single_borehole_basic.jl` - Single borehole simulation
  - `well_array_basic.jl` - Multi-well array simulation
- See the [reproducibility repository](https://github.com/cwittens/2026_DBHEs_Arrays) for advanced examples reproducing results from literature (Li et al. (2021), Hu et al. (2020), Brown et al. (2023))
