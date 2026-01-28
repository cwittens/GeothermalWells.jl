# [Methodology](@id methodology)

This page provides a brief overview of the numerical methods used in GeothermalWells.jl.
For full details, see the accompanying paper (TODO: add link when published).

## Governing Equation

GeothermalWells.jl solves the three-dimensional advection-diffusion equation for temperature ``T(x,y,z,t)``:

```math
\rho c \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) - \rho c \, v_z \frac{\partial T}{\partial z}
```

where:
- ``\rho c(x,y,z)`` is the volumetric heat capacity
- ``k(x,y,z)`` is the thermal conductivity
- ``v_z(x,y,z)`` is the vertical fluid velocity (non-zero only in the pipes)

Diffusion occurs everywhere (rock, grout, pipes, fluid), while advection acts only within the borehole pipes where fluid flows.

## Operator Splitting

The equation is solved using operator splitting, separating the problem into three subproblems:

| Operator | Method | Why |
|----------|--------|-----|
| **Vertical diffusion** | ROCK2 (stabilized explicit Runge-Kutta-Chebyshev) | Coarse vertical grid (``\Delta z \sim`` meters) allows explicit treatment |
| **Horizontal diffusion** | ADI (Alternating Direction Implicit) | Fine horizontal grid (``\Delta x, \Delta y \sim`` mm) requires implicit treatment for stability |
| **Advection** | Semi-Lagrangian | Unconditionally stable for large time steps |

This splitting is motivated by the strongly anisotropic grid: fine horizontal resolution is needed to resolve the borehole geometry, while coarser vertical resolution suffices.

## Adaptive Grid

The computational grid uses adaptive spacing:
- **Fine resolution** (``\sim 2.5`` mm) near and inside boreholes to accurately resolve the geometry
- **Coarse resolution** (up to ``\sim 10`` m) far from boreholes where gradients are weaker
- **Geometric growth** between fine and coarse regions

This allows efficient simulation of the multi-scale problem (mm-scale pipes, 100m-scale thermal influence zone, km-scale depth).

## Borehole Model

The coaxial borehole heat exchanger consists of:
- **Inner pipe**: fluid flows upward (outlet)
- **Outer annulus**: fluid flows downward (inlet)
- **Pipe walls**: steel with optional insulation on inner pipe
- **Backfill/grout**: between outer pipe and rock
- **Surrounding rock**: extends to domain boundaries

At the borehole bottom, perfect mixing is assumed: the temperature is averaged across the annulus and transferred to the inner pipe.

## GPU Acceleration

GeothermalWells.jl uses [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) for vendor-agnostic GPU kernels, enabling:
- NVIDIA GPUs via CUDA.jl
- AMD GPUs via AMDGPU.jl
- CPU execution for testing and debugging

The same code runs on all backends - just change `CPU()` to `CUDABackend()`.
