# [Methodology](@id methodology)

This page provides an overview of the numerical methods used in GeothermalWells.jl.
Also, see the accompanying paper:
[A Full Three-Dimensional GPU-Accelerated Model for Deep Borehole Heat Exchangers (DBHEs) Enabling Simulation of Well Arrays](https://pangea.stanford.edu/ERE/db/GeoConf/papers/SGW/2026/Wittenstein.pdf) or my Master Thesis for more details (TODO: add link when published).

## Governing Equation

GeothermalWells.jl solves the three-dimensional advection-diffusion equation for temperature ``T(x,y,z,t)``:

```math
\rho c \frac{\partial T}{\partial t}
= \frac{\partial}{\partial x}\!\left(k\frac{\partial T}{\partial x}\right)
  + \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)
  + \frac{\partial}{\partial z}\!\left(k\frac{\partial T}{\partial z}\right)
  - \rho c \, v_z \frac{\partial T}{\partial z}
```

where ``\rho c(x,y,z)`` is the volumetric heat capacity, ``k(x,y,z)`` is the thermal conductivity, and ``v_z(x,y,z)`` is the vertical fluid velocity (non-zero only in the borehole pipes). Diffusion occurs everywhere (rock, grout, pipes, fluid), while advection acts only within the borehole pipes where fluid flows at prescribed constant velocities.

## Adaptive Grid

The computational grid uses non-uniform Cartesian spacing to handle the multi-scale nature of the problem:

- **Fine resolution** (``\Delta x, \Delta y \sim 2.5`` mm) inside and around boreholes to resolve the centimeter-scale pipe geometry
- **Coarse resolution** (up to ``\Delta x, \Delta y \sim 10`` m) far from boreholes where temperature gradients are weaker
- **Geometric growth** between fine and coarse regions, with each successive cell a constant factor (typically 1.3) wider than the previous one
- **Uniform vertical spacing** (``\Delta z \sim 10``--``100`` m), since vertical gradients are much weaker than horizontal gradients near the borehole

For a single borehole, this typically results in about ``150 \times 150 \times 50 \approx 10^6`` grid points.

For well arrays, each borehole introduces its own set of fine grid lines. If boreholes share the same ``x``- or ``y``-coordinates (as in a regular rectangular array), they also share fine grid lines, so a ``2 \times 2`` array roughly quadruples the grid size rather than scaling worse.

## Operator Splitting

The right-hand side of the governing equation is split into three subproblems, each solved with a method tailored to its character:

```math
\begin{aligned}
  \mathcal{D}_z\, T &= \frac{1}{\rho c}\frac{\partial}{\partial z}\!\left(k\frac{\partial T}{\partial z}\right),
    &\quad&\text{(vertical diffusion)}\\[4pt]
  \mathcal{D}_{xy}\, T &= \frac{1}{\rho c}\left[\frac{\partial}{\partial x}\!\left(k\frac{\partial T}{\partial x}\right)
    + \frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right)\right],
    &\quad&\text{(horizontal diffusion)}\\[4pt]
  \mathcal{S}\, T &= -v_z\frac{\partial T}{\partial z},
    &\quad&\text{(advection)}
\end{aligned}
```

The horizontal diffusion and advection are grouped into one combined operator ``\mathcal{A}``, while the vertical diffusion forms operator ``\mathcal{B}``. These are coupled using Strang splitting: a full time step ``\Delta t`` is organized as

```math
T^{n+1} = \mathcal{A}_{\Delta t/2} \circ \mathcal{B}_{\Delta t} \circ \mathcal{A}_{\Delta t/2} \; T^n,
```

where ``\mathcal{A}_\tau`` and ``\mathcal{B}_\tau`` denote the solution operators advanced by time ``\tau``. Strang splitting is second-order accurate in ``\Delta t``, consistent with the second-order accuracy of the individual schemes.

In practice, the very first and very last half-steps of the simulation are omitted, introducing a one-time first-order error that is negligible over millions of time steps.

| Operator | Method | Rationale |
|----------|--------|-----------|
| **Vertical diffusion** ``\mathcal{B}`` | ROCK2 (stabilized explicit Runge-Kutta-Chebyshev) | Coarse vertical grid (``\Delta z \sim`` meters) allows explicit treatment |
| **Horizontal diffusion** ``\mathcal{D}_{xy}`` | ADI (Alternating Direction Implicit) | Fine horizontal grid (``\Delta x, \Delta y \sim`` mm) requires implicit treatment for stability |
| **Advection** ``\mathcal{S}`` | Semi-Lagrangian | Unconditionally stable, avoids CFL restriction from high fluid velocities |

## Vertical Diffusion: ROCK2

The vertical diffusion operator ``\mathcal{B} = \mathcal{D}_z`` is integrated using the second-order ROCK2 method of Abdulle and Medovikov, a stabilized explicit Runge-Kutta-Chebyshev scheme implemented in [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl). Standard explicit methods would require unnecessarily small time steps because their stability regions extend only a short distance along the negative real axis.

## [Horizontal Diffusion: ADI](@id adi-section)

The horizontal diffusion operator ``\mathcal{D}_{xy}`` is solved using the Peaceman-Rachford Alternating Direction Implicit (ADI) method, which splits the two-dimensional implicit problem into a sequence of one-dimensional problems. Starting from ``T^n``, one full ADI step of size ``\Delta t`` consists of:

```math
\begin{aligned}
  \left(I - \frac{\Delta t}{2}\mathcal{D}_x\right) T^* &= \left(I + \frac{\Delta t}{2}\mathcal{D}_y\right) T^n,\\[4pt]
  \left(I - \frac{\Delta t}{2}\mathcal{D}_y\right) T^{n+1} &= \left(I + \frac{\Delta t}{2}\mathcal{D}_x\right) T^*,
\end{aligned}
```

where ``\mathcal{D}_x`` and ``\mathcal{D}_y`` are the one-dimensional diffusion operators. In each half-step, one direction is treated implicitly and the other explicitly. After finite difference discretization on the Cartesian grid, each implicit direction produces a tridiagonal system that is solved in ``\mathcal{O}(N)`` operations using the Thomas algorithm. Since all tridiagonal systems along different grid lines are independent, they can be solved in parallel, which is particularly beneficial on GPU architectures.

The structured Cartesian grid is what guarantees the tridiagonal structure: every grid line in the ``x``-direction contains the same set of ``x``-coordinates regardless of the ``y``- and ``z``-position. On an unstructured grid, the implicit systems would generally be sparse but not tridiagonal.

## Advection: Semi-Lagrangian Method

The advection operator ``\mathcal{S}`` transports heat along the borehole pipes at prescribed constant fluid velocities: ``v_\text{outer}`` (downward in the outer annulus) and ``v_\text{inner}`` (upward in the inner pipe). A standard explicit upwind scheme would impose a CFL condition ``\Delta t \leq \Delta z / |v_z|``, restricting time steps to a few seconds.

The semi-Lagrangian method avoids this by tracing characteristic lines backward in time. The updated temperature at grid point ``z_i`` is

```math
T^{n+1}(z_i) = T^n(z_i - v_z \Delta t),
```

where the departure point ``z_i - v_z \Delta t`` generally does not coincide with a grid point. The temperature there is obtained by linear interpolation between the two bracketing grid points:

```math
T^n(z_i - v_z \Delta t) = (1 - \alpha)\, T^n(z_k) + \alpha\, T^n(z_{k+1}).
```

This method is unconditionally stable with respect to the advection velocity, allowing arbitrarily large time steps. At the borehole bottom (depth ``h``), where fluid transitions from the outer annulus to the inner pipe, perfect mixing is assumed: the temperature is averaged across the annulus cross-section and transferred to the inner pipe.

## Putting It Together

Operator ``\mathcal{A}`` combines horizontal diffusion and advection internally using two ADI half-steps with semi-Lagrangian advection interleaved. Denoting the advection operator advanced by ``\tau`` as ``\mathcal{S}_\tau``, a single application of ``\mathcal{A}_\tau`` proceeds as:

**First ADI half-step** (``y``-explicit, ``x``-implicit):
```math
\begin{aligned}
  T^{(1)} &= \bigl(I + \tfrac{\tau}{2}\,\mathcal{D}_y\bigr)\, T^{n},\\
  T^{(2)} &= \mathcal{S}_{\tau/2}\, T^{(1)},\\
  \bigl(I - \tfrac{\tau}{2}\,\mathcal{D}_x\bigr)\, T^{(3)} &= T^{(2)}.
\end{aligned}
```

**Second ADI half-step** (``x``-explicit, ``y``-implicit):
```math
\begin{aligned}
  T^{(4)} &= \bigl(I + \tfrac{\tau}{2}\,\mathcal{D}_x\bigr)\, T^{(3)},\\
  T^{(5)} &= \mathcal{S}_{\tau/2}\, T^{(4)},\\
  \bigl(I - \tfrac{\tau}{2}\,\mathcal{D}_y\bigr)\, \widetilde{T} &= T^{(5)},
\end{aligned}
```

where ``\widetilde{T} = \mathcal{A}_\tau\, T^n``. The advection is placed after the explicit diffusion update and before the implicit Thomas solve. This ordering was determined empirically; placing the advection elsewhere led to numerical instabilities in all tested configurations.

Typical full time steps are on the order of ``\Delta t \approx 50``--``150`` s, several orders of magnitude larger than what a naive explicit treatment of the horizontal diffusion or advection would allow.

## Borehole Model

The coaxial borehole heat exchanger consists of:
- **Inner pipe**: fluid flows upward (toward outlet)
- **Outer annulus**: fluid flows downward (from inlet)
- **Pipe walls**: steel, with optional insulation on the inner pipe to reduce thermal short-circuiting
- **Backfill/grout**: between outer pipe and rock
- **Surrounding rock**: extends to domain boundaries

At the borehole bottom, perfect mixing is assumed: the temperature is averaged across the annulus and transferred to the inner pipe.

All spatial derivatives in the diffusion operators are approximated using standard second-order finite differences with variable coefficients on the non-uniform grid.

## GPU Acceleration

GeothermalWells.jl uses [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) for vendor-agnostic GPU kernels, enabling:
- NVIDIA GPUs via CUDA.jl
- AMD GPUs via AMDGPU.jl
- CPU execution for testing and debugging

The same code runs on all backends - just change `CPU()` to `CUDABackend()`.
