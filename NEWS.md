# Changelog

GeothermalWells.jl follows [semantic versioning](https://semver.org/).

## Unreleased

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
