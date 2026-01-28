# [API Reference](@id api-reference)

## Material Properties

```@docs
HomogenousMaterialProperties
StratifiedMaterialProperties
```

## Borehole Geometry

```@docs
Borehole
```

## Inlet Models

```@docs
ConstantInlet
HeatExchangerInlet
CustomInlet
```

## Grid Generation

```@docs
create_adaptive_grid_1d
create_uniform_gridz_with_borehole_depths
compute_domain
```

## Initial Conditions

```@docs
initial_condition_thermal_gradient
```

## Simulation

```@docs
create_cache
rhs_diffusion_z!
eigen_estimator
get_simulation_callback
```

## Utilities

```@docs
pkg_dir
examples_dir
data_dir
plot_grid
```

## Validation Data

Functions for loading reference data from published studies:

```@docs
data_li
data_hu
data_brown_single_well_b
data_brown_single_well_c
```
