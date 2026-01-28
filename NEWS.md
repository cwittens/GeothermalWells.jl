# Changelog

GeothermalWells.jl follows [semantic versioning](https://semver.org/).

## Changes in v0.1.1

- renamed `get_callback` to `get_simulation_callback` to make clear its needed and not optional ([#4]).
- added Documentation ([#4]). 
- Released as [official Julia package](https://github.com/JuliaRegistries/General/pull/146869)
- Bug fix with `print_every_n` parameter in callback

## v0.1.0

Initial release with core functionality:
- Coaxial borehole heat exchanger simulation
- Adaptive grid generation
- GPU support via KernelAbstractions.jl
- Integration with OrdinaryDiffEq.jl solvers
