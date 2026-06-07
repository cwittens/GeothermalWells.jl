# Changelog

GeothermalWells.jl follows [semantic versioning](https://semver.org/).

## Unreleased

## v0.3.1

### Added
- Added Brown et al. validation data loaders for Figure 6 array cases: `data_brown_array(spacing, i)`  

### Changed
- `data_brown_single_well_c` had the wrong name and is now `data_brown_single_well_d` (but `data_brown_single_well_c` is still there for legacy-compatible)

## v0.3.0

### Added
- `create_cache` now accepts a `precompute_materials` keyword argument controlling how material properties are looked up inside the diffusion kernels. Accepted values: `:auto` (default; picks `:on_the_fly` for a single borehole and `:precomputed` for well arrays), `true`/`:precomputed`, or `false`/`:on_the_fly`/`:onthefly`. Precomputed material arrays are significantly faster for multi-borehole simulations; on-the-fly evaluation is slightly faster for a single well.
- New `AbstractMaterialAccessor` interface with `PrecomputedMaterialAccessor` and `OnTheFlyMaterialAccessor` structs used internally by the diffusion kernels.
- `Float_used` is now optional and defaults to `Float64`:
  - `create_adaptive_grid_1d(...; Float_used=Float64, ...)` (was a required keyword argument).
  - `initial_condition_thermal_gradient(backend, gridx, gridy, gridz; T_surface, gradient, Float_used=Float64)` — new 4-positional form. The legacy 5-positional form `initial_condition_thermal_gradient(backend, Float_used, gridx, gridy, gridz; ...)` is kept for backward compatibility.
  - `Borehole(xc, yc, h, r_inner, t_inner, r_outer, t_outer, r_backfill, ṁ, insulation_depth)` — outer constructor defaulting to `Float64`. `Borehole{T}(...)` still works for any `T<:Real`.

## v0.2.1

### Changed
- updated restart internals to address edge case

## v0.2.0

### Added
- Checkpointing and restart via `prepare_restart`, `reload_snapshots!`, and new `checkpoint_dir`/`checkpoint_id`/`checkpoint_every_n` keyword arguments on `get_simulation_callback`
- Dirichlet-like BC for vertical diffusion at the surface

### Changed
- **Breaking:** Operator splitting changed from Lie splitting to Strang splitting for improved accuracy
- **Breaking:** Removed `write_to_jld` and `data_folder_dir` keyword arguments from `get_simulation_callback`
- Advection now precomputes turnaround mean temperature, improving efficiency for well array simulations
- Improved documentation on methodology

## v0.1.3

### Added
- Logo and favicon to documentation

### Changed
- Added validation error when `write_to_jld=true` but `data_folder_dir` is empty

## v0.1.2

### Added
- Documentation ([#4])

### Changed
- Renamed `get_callback` to `get_simulation_callback` to clarify that it is required and not optional ([#4])

---

## v0.1.1

### Added
- Released as an official Julia package in the General registry

### Fixed
- Bug in `print_every_n` parameter handling in the callback

---

## v0.1.0

### Added
- Initial release with core functionality
- Coaxial borehole heat exchanger simulation
- Adaptive grid generation
- GPU support via KernelAbstractions.jl
- Integration with OrdinaryDiffEq.jl solvers
