module GeothermalWells

using Adapt: adapt
# using GPUArraysCore: @allowscalar
using KernelAbstractions: @kernel, @index, @Const, @uniform, @private, @atomic, CPU
import KernelAbstractions: zeros
using DiffEqCallbacks: SavedValues, SavingCallback, CallbackSet
using MuladdMacro: @muladd
using DiffEqBase: DiscreteCallback
using DelimitedFiles: readdlm
using JLD2: @save
using Plots: plot, plot!, vline!, hline!, annotate!, scatter, current, text

include("physics.jl")
include("cache.jl")
include("grids.jl")
include("solver.jl")
include("util.jl")
# include some plotting functions

# Exports from physics.jl
export AbstractMaterialProperties,
    AbstractBorehole,
    Borehole,
    StratifiedMaterialProperties,
    HomogenousMaterialProperties,
    ConstantInlet,
    HeatExchangerInlet,
    CustomInlet,
    initial_condition_thermal_gradient

# Exports from cache.jl
export create_cache,
    get_callback

# Exports from grids.jl
export compute_domain,
    create_uniform_gridz_with_borehole_depths,
    create_adaptive_grid_1d

# Exports from solver.jl
export rhs_diffusion_z!

# Exports from util.jl
export pkg_dir,
    examples_dir,
    data_dir,
    data_brown_single_well_b,
    data_brown_single_well_c,
    data_hu,
    data_li,
    plot_grid
end
