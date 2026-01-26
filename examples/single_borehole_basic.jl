# Single deep borehole heat exchanger simulation
# Set up inspired from: Hu et al. (2020) - Numerical modeling of a coaxial borehole heat exchanger
# https://doi.org/10.1016/j.renene.2019.09.141

using GeothermalWells
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
using KernelAbstractions: CPU

# Choose backend: CPU() for testing, or CUDABackend()/ROCBackend() for GPU
backend = CPU()
Float_used = Float64

# =============================================================================
# Material properties (homogeneous rock)
# =============================================================================
materials = HomogenousMaterialProperties{Float_used}(
    2.88,       # k_rock [W/(m·K)]
    2.17e6,     # rho_c_rock [J/(m³·K)]
    0.6,        # k_water [W/(m·K)]
    4.18e6,     # rho_c_water [J/(m³·K)]
    44.5,       # k_steel [W/(m·K)]
    3.73e6,     # rho_c_steel [J/(m³·K)]
    0.26,       # k_insulating [W/(m·K)]
    1.96e6,     # rho_c_insulating [J/(m³·K)]
    1.0,        # k_backfill [W/(m·K)] (not used - no backfill)
    1.0         # rho_c_backfill [J/(m³·K)] (not used)
)

# =============================================================================
# Borehole geometry
# =============================================================================
# Deep coaxial borehole heat exchanger (350m depth) shallow to make it faster to simulate on CPU
# Inner pipe is insulated down to 1000m to reduce thermal short-circuiting
borehole = Borehole{Float_used}(
    0.0,             # xc [m]
    0.0,             # yc [m]
    350.0,          # h - borehole depth [m]
    0.0381,          # r_inner - inner pipe radius [m]
    0.01,            # t_inner - inner pipe wall thickness [m]
    0.0889,          # r_outer - outer pipe inner radius [m]
    0.01,            # t_outer - outer pipe wall thickness [m]
    0.0989,          # r_backfill - outer radius (r_outer + t_outer) [m]
    10.0,            # ṁ - mass flow rate [kg/s]
    100.0           # insulation_depth [m]
)

boreholes = (borehole,)

# =============================================================================
# Grid setup
# =============================================================================
# Domain boundaries
xmin, xmax = -100.0, 100.0
ymin, ymax = -100.0, 100.0
zmin, zmax = 0.0, 370.0

# Grid parameters
dx_fine = 0.0025      # fine spacing near borehole [m]
growth_factor = 1.3   # geometric growth rate
dx_max = 10.0         # maximum spacing far from borehole [m]
dz = 10.0            # vertical spacing [m]

# Create adaptive grids (fine near borehole, coarse far away)
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

gridz = create_uniform_gridz_with_borehole_depths(zmin=zmin, zmax=zmax, dz=dz, boreholes=boreholes, backend=backend)

println("Grid size: $(length(gridx)) x $(length(gridy)) x $(length(gridz))")

# =============================================================================
# Initial condition
# =============================================================================
# Linear thermal gradient: T(z) = T_surface + gradient * z
T0 = initial_condition_thermal_gradient(
    backend, Float_used, gridx, gridy, gridz;
    T_surface=2.29,    # surface temperature [°C]
    gradient=0.35     # thermal gradient [K/m] (10x higher than typical to speed up testing on CPU)
);


# =============================================================================
# Inlet model
# =============================================================================
# Constant inlet temperature of 20°C
inlet_model = ConstantInlet{Float_used}(20.0)

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
tspan = (0.0, 800.0) 

prob = ODEProblem(rhs_diffusion_z!, T0, tspan, cache)

# Save solution at regular intervals
n_saves = 11  # save initial + 10 more times
saveat = range(tspan..., n_saves)
callback, saved_values = get_callback(
    saveat=saveat,
    print_every_n=1,
    write_to_jld=false
)

# Time step and solver
Δt = 80.0  # [s]

println("Simulating with Δt = $(Δt)s")

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