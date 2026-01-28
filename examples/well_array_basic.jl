# Well array simulation (2x2 array of deep borehole heat exchangers)
# This example demonstrates how to simulate multiple boreholes arranged in a grid
# For GPU acceleration, change backend to CUDABackend() or ROCBackend()

using GeothermalWells
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
using KernelAbstractions: CPU

# Choose backend: CPU() for testing, or CUDABackend()/ROCBackend() for GPU
backend = CPU()
Float_used = Float64

# =============================================================================
# Borehole array configuration
# =============================================================================
# Define array layout by specifying x and y coordinates for borehole centers
# 2x2 array with 30m spacing (small array for CPU testing)
XC = [-15.0, 15.0]
YC = [-15.0, 15.0]

# Alternative configurations (uncomment to use):
# 3x3 array with 30m spacing
# XC = [-30.0, 0.0, 30.0]
# YC = [-30.0, 0.0, 30.0]

# 5x5 array with 30m spacing (use GPU for larger arrays)
# XC = [-60.0, -30.0, 0.0, 30.0, 60.0]
# YC = [-60.0, -30.0, 0.0, 30.0, 60.0]

println("Configuring $(length(XC))x$(length(YC)) borehole array")

# =============================================================================
# Material properties (homogeneous rock)
# =============================================================================
materials = HomogenousMaterialProperties{Float_used}(
    2.5,        # k_rock [W/(m·K)]
    2.0e6,      # rho_c_rock [J/(m³·K)]
    0.6,        # k_water [W/(m·K)]
    4.18e6,     # rho_c_water [J/(m³·K)]
    44.5,       # k_steel [W/(m·K)]
    3.73e6,     # rho_c_steel [J/(m³·K)]
    0.26,       # k_insulating [W/(m·K)]
    1.96e6,     # rho_c_insulating [J/(m³·K)]
    1.5,        # k_backfill [W/(m·K)]
    1.76e6      # rho_c_backfill [J/(m³·K)]
)

# =============================================================================
# Borehole geometry
# =============================================================================
# Create array of boreholes at specified coordinates
# All boreholes have identical geometry (can be customized if needed)
# Using shallow depth (350m) for faster CPU testing
boreholes = tuple(
    (Borehole{Float_used}(
        xc,                      # xc [m]
        yc,                      # yc [m]
        350.0,                   # h - borehole depth [m] (shallow for CPU testing)
        0.0381,                  # r_inner - inner pipe radius [m]
        0.01,                    # t_inner - inner pipe wall thickness [m]
        0.0889,                  # r_outer - outer pipe inner radius [m]
        0.01,                    # t_outer - outer pipe wall thickness [m]
        0.0989,                  # r_backfill - outer radius [m]
        10.0,                    # ṁ - mass flow rate [kg/s]
        100.0                    # insulation_depth [m]
    ) for xc in XC, yc in YC)...
)

println("Created $(length(boreholes)) boreholes")

# =============================================================================
# Grid setup
# =============================================================================
# Domain boundaries (adjusted to encompass entire array with buffer)
xmin, xmax = -80.0, 80.0
ymin, ymax = -80.0, 80.0
zmin, zmax = 0.0, 370.0

# Grid parameters
dx_fine = 0.0025      # fine spacing near each borehole [m]
growth_factor = 1.3   # geometric growth rate
dx_max = 10.0         # maximum spacing far from boreholes [m]
dz = 10.0             # vertical spacing [m]

# Create adaptive grids (fine near all boreholes, coarse elsewhere)
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

# =============================================================================
# Initial condition
# =============================================================================
# Linear thermal gradient: T(z) = T_surface + gradient * z
T0 = initial_condition_thermal_gradient(
    backend, Float_used, gridx, gridy, gridz;
    T_surface=10.0,     # surface temperature [°C]
    gradient=0.35       # thermal gradient [K/m] (10x higher than typical to speed 
                        # up testing on CPU, by needing smaller zmax)
)

# =============================================================================
# Inlet model
# =============================================================================
# Heat exchanger inlet: all boreholes extract heat with same ΔT
# T_inlet = T_outlet - ΔT
ΔT = 5.0  # temperature difference [K]
inlet_model = HeatExchangerInlet{Float_used}(ΔT)

# Alternative: constant inlet temperature for all boreholes
# inlet_model = ConstantInlet{Float_used}(20.0)

# =============================================================================
# Create simulation cache
# =============================================================================
cache = create_cache(
    backend=backend,
    gridx=gridx,
    gridy=gridy,
    gridz=gridz,
    materials=materials,
    boreholes=boreholes,
    inlet_model=inlet_model
)

# =============================================================================
# Time integration
# =============================================================================
tspan = (0.0, 800.0)  # short simulation for testing [s]

prob = ODEProblem(rhs_diffusion_z!, T0, tspan, cache)

# Save solution at regular intervals
n_saves = 5
saveat = range(tspan..., n_saves)
callback, saved_values = get_simulation_callback(
    saveat=saveat,
    print_every_n=1
)

# Time step and solver
Δt = 80.0  # [s]

println("Simulating $(length(XC))x$(length(YC)) array with Δt = $(Δt)s")

t_elapsed = @elapsed solve(
    prob,
    ROCK2(max_stages=100, eigen_est=eigen_estimator),
    save_everystep=false,
    callback=callback,
    adaptive=false,
    dt=Δt,
    maxiters=Int(1e10)
)

println("Simulation completed in $(round(t_elapsed, digits=2))s")
