# Changelog

GeothermalWells.jl follows [semantic versioning](https://semver.org/).

## Unreleased

## v0.2.0

### Added
- Fault-tolerant checkpointing and restart via `prepare_restart`, `reload_snapshots!`, and new `checkpoint_dir`/`checkpoint_id`/`checkpoint_every_n` keyword arguments on `get_simulation_callback`
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
