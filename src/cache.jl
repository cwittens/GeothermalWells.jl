"""
    AbstractMaterialAccessor

Small kernel-facing abstraction for material-property lookup.

`PrecomputedMaterialAccessor` stores precomputed thermal-conductivity and
volumetric-heat-capacity arrays. `OnTheFlyMaterialAccessor` stores the borehole
geometry and material model and evaluates the existing geometry-dependent
material functions inside the kernel.
"""
abstract type AbstractMaterialAccessor end

struct PrecomputedMaterialAccessor{KArr,RhoCArr} <: AbstractMaterialAccessor
    thermal_conductivity::KArr
    volumetric_heat_capacity::RhoCArr
end

struct OnTheFlyMaterialAccessor{B,M} <: AbstractMaterialAccessor
    boreholes::B
    materials::M
end

@adapt_structure PrecomputedMaterialAccessor
@adapt_structure OnTheFlyMaterialAccessor

@inline function lookup_thermal_conductivity(mat::PrecomputedMaterialAccessor, i, j, k, gridx, gridy, gridz)
    return mat.thermal_conductivity[k, j, i]
end

@inline function lookup_volumetric_heat_capacity(mat::PrecomputedMaterialAccessor, i, j, k, gridx, gridy, gridz)
    return mat.volumetric_heat_capacity[k, j, i]
end

@inline function lookup_thermal_conductivity(mat::OnTheFlyMaterialAccessor, i, j, k, gridx, gridy, gridz)
    return get_thermal_conductivity(gridx[i], gridy[j], gridz[k], mat.boreholes, mat.materials)
end

@inline function lookup_volumetric_heat_capacity(mat::OnTheFlyMaterialAccessor, i, j, k, gridx, gridy, gridz)
    return get_volumetric_heat_capacity(gridx[i], gridy[j], gridz[k], mat.boreholes, mat.materials)
end

@kernel function precompute_materials_kernel!(k_arr, rho_c_arr, @Const(gridx), @Const(gridy), @Const(gridz), boreholes, materials)
    k, j, i = @index(Global, NTuple)
    x, y, z = gridx[i], gridy[j], gridz[k]
    k_arr[k, j, i] = get_thermal_conductivity(x, y, z, boreholes, materials)
    rho_c_arr[k, j, i] = get_volumetric_heat_capacity(x, y, z, boreholes, materials)
end


function _material_mode(precompute_materials, N_bh)
    if precompute_materials === :auto
        return N_bh > 1 ? :precomputed : :on_the_fly

    elseif precompute_materials == true || precompute_materials == :precomputed
        return :precomputed

    elseif precompute_materials == false ||
           precompute_materials == :on_the_fly ||
           precompute_materials == :onthefly
        return :on_the_fly

    else
        throw(ArgumentError(
            "precompute_materials must be :auto, true/false, :precomputed, or :on_the_fly"
        ))
    end
end


"""
    create_advection_index_lists(backend, gridx, gridy, gridz, boreholes)

Create index lists for advection in inner and outer pipes.

Returns tuple with index lists for efficient advection kernel dispatch.
"""
function create_advection_index_lists(backend, gridx, gridy, gridz, boreholes)

    # calculate index map (grouped by borehole, in borehole order)
    Idx_list_Inner = Vector{Tuple{Int,Int,Int}}()
    Idx_list_Outer = Vector{Tuple{Int,Int,Int}}()
    gridx_cpu = adapt(CPU(), gridx)
    gridy_cpu = adapt(CPU(), gridy)

    for (bh_idx, bh) in enumerate(boreholes)
        for (i, x) in enumerate(gridx_cpu), (j, y) in enumerate(gridy_cpu)
            r_sq = (x - bh.xc)^2 + (y - bh.yc)^2
            if r_sq < bh.r_inner^2
                push!(Idx_list_Inner, (i, j, bh_idx))
            elseif (bh.r_inner + bh.t_inner)^2 <= r_sq < bh.r_outer^2
                push!(Idx_list_Outer, (i, j, bh_idx))
            end
        end
    end

    # number of inner/outer cells of each borehole and the offsets of each
    # borehole's (contiguous) block within the index lists
    N_bh = length(boreholes)
    count_inner_per_bh_cpu = Base.zeros(Int, N_bh)
    for (_, _, bh_idx) in Idx_list_Inner
        count_inner_per_bh_cpu[bh_idx] += 1
    end
    count_outer_per_bh_cpu = Base.zeros(Int, N_bh)
    for (_, _, bh_idx) in Idx_list_Outer
        count_outer_per_bh_cpu[bh_idx] += 1
    end
    inner_offset_per_bh = adapt(backend, cumsum(count_inner_per_bh_cpu) .- count_inner_per_bh_cpu)
    outer_offset_per_bh = adapt(backend, cumsum(count_outer_per_bh_cpu) .- count_outer_per_bh_cpu)
    count_inner_per_bh = adapt(backend, count_inner_per_bh_cpu)
    count_outer_per_bh = adapt(backend, count_outer_per_bh_cpu)

    countxy_inner = length(Idx_list_Inner)
    countxy_outer = length(Idx_list_Outer)
    Idx_list = vcat(Idx_list_Inner, Idx_list_Outer)
    Idx_list = adapt(backend, Idx_list)
    Idx_list_Outer = adapt(backend, Idx_list_Outer)
    Idx_list_Inner = adapt(backend, Idx_list_Inner)

    # just take the maximum h for countz
    # eg have a little extra space (if some boreholes are shorter)
    # but therefore simpler indexing...
    max_h = maximum(bh.h for bh in boreholes)
    countz = sum(gridz .<= max_h)

    # the fused advection kernel stages one (i,j) column of ϕ in shared memory and
    # needs a workgroup spanning all countz levels; fall back to the two-kernel
    # u_tmp version when countz exceeds the maximum workgroup size
    if countz <= MAX_FUSED_ADVECTION_COUNTZ
        u_tmp = nothing
    else
        # layout (countz, countxy): consecutive threads (fastest in k) access consecutive memory
        u_tmp = zeros(backend, eltype(gridx), countz, countxy_inner + countxy_outer)
    end

    return (Idx_list_Inner, Idx_list_Outer, Idx_list,
        count_inner_per_bh, count_outer_per_bh, inner_offset_per_bh, outer_offset_per_bh,
        countxy_inner, countxy_outer, countz, u_tmp)
end

"""
    create_cache(; backend, gridx, gridy, gridz, materials, boreholes, inlet_model,
                   precompute_materials=:auto, precompute_thomas=true)

Create simulation cache.

By default, `precompute_materials=:auto` chooses the material evaluation mode based
on the number of boreholes:

- one borehole: use `:on_the_fly`
- multiple boreholes: use `:precomputed`

This reflects the observed performance behavior: on-the-fly material lookup is slightly
faster for a single well, while precomputed material arrays are significantly faster
for well arrays.

You can override the automatic choice manually:

- `precompute_materials=true` or `:precomputed`
- `precompute_materials=false` or `:on_the_fly`

With `precompute_thomas=true` (the default), the tridiagonal factorizations of the ADI
implicit solves are precomputed once (they depend only on the constant time step, the
grid, and the time-independent material distribution) and each implicit solve becomes a
lean two-sweep kernel. This costs six additional arrays of the size of the temperature
field but speeds up the dominant ADI solves considerably. Set `precompute_thomas=false`
to rebuild the systems in every solve (lower memory footprint, original behavior).

Returns named tuple containing grids, materials, index lists, outlet temperature arrays,
eigenvalue estimates, material accessor, and precomputed `Val` types for kernel dispatch.
"""
function create_cache(; backend, gridx, gridy, gridz, materials, boreholes, inlet_model, precompute_materials=:auto, precompute_thomas=true)

    Nx, Ny, Nz = length(gridx), length(gridy), length(gridz)
    N_bh = length(boreholes)
    material_mode = _material_mode(precompute_materials, N_bh)

    Idx_list_Inner, Idx_list_Outer, Idx_list,
    count_inner_per_bh, count_outer_per_bh, inner_offset_per_bh, outer_offset_per_bh,
    countxy_inner, countxy_outer, countz, u_tmp = create_advection_index_lists(backend, gridx, gridy, gridz, boreholes)

    T_outlet = zeros(backend, eltype(gridx), N_bh)

    T_turnaround_mean = zeros(backend, eltype(gridx), countz, N_bh)

    eigen_estimate = eigen_estimator_pre_calculation(gridz, materials)

    # pre compute Val types for kernels
    ValNx = Val(Nx)
    ValNy = Val(Ny)
    ValNz = Val(Nz)
    ValTrue = Val(true)
    ValFalse = Val(false)
    Val_in_x = Val(:x)
    Val_in_y = Val(:y)
    Val_in_z = Val(:z)

    gridz_cpu = adapt(CPU(), gridz)

    for bh in boreholes
        if !(bh.h in gridz_cpu)
            throw(ArgumentError("Borehole depth h=$(bh.h) must be present in gridz. Use create_uniform_gridz_with_borehole_depths or include each borehole depth explicitly."))
        end
    end

    gridx = adapt(backend, gridx)
    gridy = adapt(backend, gridy)
    gridz = adapt(backend, gridz)

    if material_mode === :precomputed
        Thermal_Conductivity = zeros(backend, eltype(gridx), Nz, Ny, Nx)
        Volumetric_Heat_Capacity = zeros(backend, eltype(gridx), Nz, Ny, Nx)
        precompute_materials_kernel!(backend)(Thermal_Conductivity, Volumetric_Heat_Capacity, gridx, gridy, gridz, boreholes, materials, ndrange=(Nz, Ny, Nx))
        material_accessor = PrecomputedMaterialAccessor(Thermal_Conductivity, Volumetric_Heat_Capacity)
    elseif material_mode === :on_the_fly
        material_accessor = OnTheFlyMaterialAccessor(boreholes, materials)
    end

    # per-borehole scalar properties as flat arrays for the advection kernels
    # (dynamically indexing into the boreholes tuple inside a kernel generates very
    # inefficient GPU code that gets worse with the number of boreholes)
    Float_used = eltype(gridx)
    bh_h = adapt(backend, Float_used[bh.h for bh in boreholes])
    bh_v_inner = adapt(backend, Float_used[bh.v_inner for bh in boreholes])
    bh_v_outer = adapt(backend, Float_used[bh.v_outer for bh in boreholes])

    # storage for the precomputed Thomas factorization of the ADI implicit solves;
    # filled lazily by update_thomas_factors! once the time step is known
    if precompute_thomas
        thomas_factors = (;
            dt=Ref(NaN), # Float64 so the comparison with the integrator's dt is exact
            W_x=zeros(backend, Float_used, Nz, Ny, Nx),
            invD_x=zeros(backend, Float_used, Nz, Ny, Nx),
            Uinv_x=zeros(backend, Float_used, Nz, Ny, Nx),
            W_y=zeros(backend, Float_used, Nz, Ny, Nx),
            invD_y=zeros(backend, Float_used, Nz, Ny, Nx),
            Uinv_y=zeros(backend, Float_used, Nz, Ny, Nx),
        )
    else
        thomas_factors = nothing
    end

    cache = (;
        backend,
        material_accessor,
        gridx,
        gridy,
        gridz,
        Nx,
        Ny,
        Nz,
        N_bh,
        materials,
        boreholes,
        inlet_model,
        T_outlet,
        T_turnaround_mean,
        u_tmp,
        Idx_list_Inner,
        Idx_list_Outer,
        Idx_list,
        count_inner_per_bh,
        count_outer_per_bh,
        inner_offset_per_bh,
        outer_offset_per_bh,
        countxy_inner,
        countxy_outer,
        countz,
        Val_countz=Val(countz),
        bh_h,
        bh_v_inner,
        bh_v_outer,
        thomas_factors,
        eigen_estimate,
        ValNx,
        ValNy,
        ValNz,
        ValTrue,
        ValFalse,
        Val_in_x,
        Val_in_y,
        Val_in_z,
    )

    return cache
end





# =============================================================================
# Checkpoint and snapshot file helpers
# =============================================================================
 
"""
    _checkpoint_path(checkpoint_dir, checkpoint_id) -> String
 
Return path for the restart checkpoint file: `checkpoint_dir/checkpoint_{id}.jld2`.
"""
_checkpoint_path(checkpoint_dir, checkpoint_id) = joinpath(checkpoint_dir, "checkpoint_$(checkpoint_id).jld2")
 
"""
    _snapshot_path(checkpoint_dir, checkpoint_id, n) -> String
 
Return path for snapshot number `n`: `checkpoint_dir/snapshot_{id}_{NNNN}.jld2`.
"""
_snapshot_path(checkpoint_dir, checkpoint_id, n) = joinpath(checkpoint_dir, "snapshot_$(checkpoint_id)_$(lpad(n, 4, '0')).jld2")

function _is_snapshot_file(filename, checkpoint_id)
    prefix = "snapshot_$(checkpoint_id)_"
    suffix = ".jld2"

    startswith(filename, prefix) || return false
    endswith(filename, suffix) || return false

    number_part_with_suffix = filename[nextind(filename, lastindex(prefix)):end]
    number_part = chop(number_part_with_suffix; tail=length(suffix))
    return !isempty(number_part) && all(isdigit, number_part)
end
 
"""
    _load_existing_snapshots(checkpoint_dir, checkpoint_id, Float_used_to_save)
 
Scan `checkpoint_dir` for existing snapshot files matching `checkpoint_id`, load them
in chronological order, and return `(times, arrays, count)`.
"""
function _load_existing_snapshots(checkpoint_dir, checkpoint_id, Float_used_to_save)
    times = Float64[]
    arrays = Array{Float_used_to_save, 3}[]
 
    if !isdir(checkpoint_dir)
        return times, arrays, 0
    end
 
    files = filter(f -> _is_snapshot_file(f, checkpoint_id), readdir(checkpoint_dir))
    sort!(files)
 
    for f in files
        path = joinpath(checkpoint_dir, f)
        @load path u_save t_save
        push!(times, Float64(t_save))
        push!(arrays, Float_used_to_save.(u_save))
    end
 
    return times, arrays, length(files)
end
 
"""
    _clean_and_count_snapshots(checkpoint_dir, checkpoint_id) -> Int
 
Determine the correct snapshot counter by checking consistency with the checkpoint state.
 
**No checkpoint file exists:** All existing snapshot files for this ID are stale (either from
a previous simulation with different parameters, or from a crashed run with no checkpoint).
Since the simulation will start from `t=0` and regenerate all saveat times, the old files
are deleted and the counter starts at 0.
 
**Checkpoint file exists at `t_checkpoint`:** Snapshots with `t_save <= t_checkpoint` are
valid (produced before the checkpoint). Snapshots with `t_save > t_checkpoint` are stale
(produced after the checkpoint in a run that later crashed, so the checkpoint doesn't
reflect them). Stale files are deleted and the counter is set to the number of valid files.
"""
function _clean_and_count_snapshots(checkpoint_dir, checkpoint_id)
    if !isdir(checkpoint_dir)
        return 0
    end
 
    files = filter(f -> _is_snapshot_file(f, checkpoint_id), readdir(checkpoint_dir))
    sort!(files)
 
    if isempty(files)
        return 0
    end
 
    # Check for checkpoint
    cp_path = _checkpoint_path(checkpoint_dir, checkpoint_id)
    if isfile(cp_path)
        @load cp_path t_checkpoint
 
        # Keep snapshots with t_save <= t_checkpoint, delete the rest
        valid_count = 0
        for f in files
            path = joinpath(checkpoint_dir, f)
            @load path t_save
            if t_save <= t_checkpoint
                valid_count += 1
            else
                rm(path)
                println("  Removed stale snapshot $(f) (t=$(t_save) > t_checkpoint=$(t_checkpoint))")
            end
        end
        return valid_count
    else
        # No checkpoint: starting from t=0, all existing snapshots are stale
        for f in files
            rm(joinpath(checkpoint_dir, f))
        end
        println("  Removed $(length(files)) stale snapshot files (no checkpoint found, starting fresh).")
        return 0
    end
end
 
 
# =============================================================================
# Callbacks
# =============================================================================
 
"""
    get_simulation_callback(; saveat, print_every_n=1000,
                              checkpoint_dir="", checkpoint_id="latest",
                              checkpoint_every_n=0, Float_used_to_save=Float32)
 
Create the required callback set for the simulation.
 
!!! warning "Required"
    This callback is essential for the simulation. It performs the ADI (Alternating Direction
    Implicit) method for horizontal diffusion and the semi-Lagrangian advection. Without this
    callback, only vertical diffusion (handled by ROCK2) is computed.
 
The callback combines up to four components:
1. **ADI + Advection**: Horizontal diffusion and fluid advection (runs every timestep)
2. **Checkpointing** (optional): Periodically saves `(u, t)` to disk for fault-tolerant restarts
3. **Progress printing**: Prints simulation progress every `print_every_n` steps
4. **Solution saving**: Saves temperature field at times specified by `saveat`. When
   `checkpoint_dir` is provided, snapshots are also written to disk so that data persists
   across crash/restart cycles.
 
The checkpoint callback is placed directly after the ADI + advection callback in the
`CallbackSet` ordering. This ensures that a checkpoint always represents a fully consistent
state (both ROCK2 vertical diffusion and ADI horizontal diffusion + advection have been
applied), so restarting from a checkpoint produces identical results to an uninterrupted run.
 
!!! note "Post-solve snapshot assembly"
    The ODE solver clears `saved_values` at the start of `solve()`, so it only contains
    snapshots from the current run. To get the full history across all crash/restart cycles,
    call [`reload_snapshots!`](@ref) after the solve completes.
 
# Arguments
- `saveat`: Times at which to save the solution (e.g., `range(0, 3600, 10)` or `[0.0, 3600.0]`)
- `print_every_n=1000`: Print progress every N accepted timesteps
- `checkpoint_dir=""`: Directory for checkpoint and snapshot files. Empty string disables
    both checkpointing and persistent snapshots (solutions are only kept in memory).
- `checkpoint_id="latest"`: Unique identifier for checkpoint/snapshot files (useful when
    multiple simulations share the same directory). Files are named
    `checkpoint_{id}.jld2` (restart) and `snapshot_{id}_NNNN.jld2` (data).
- `checkpoint_every_n=0`: Save a restart checkpoint every N accepted timesteps. Set to 0
    to disable restart checkpointing (snapshots at `saveat` times are still written if
    `checkpoint_dir` is provided).
- `Float_used_to_save=Float32`: Floating point type for saved solution snapshots. Restart
    checkpoints always use full `Float64` precision.
 
# Returns
- `callback`: Combined `CallbackSet` to pass to `solve(..., callback=callback)`
- `saved_values`: `SavedValues` object. After `solve`, contains only snapshots from the
  current run. Call [`reload_snapshots!`](@ref) to populate with the full history from disk.
 
# Example
```julia
callback, saved_values = get_simulation_callback(
    saveat=saveat,
    print_every_n=100_000,
    checkpoint_dir="output/",
    checkpoint_id="my_simulation",
    checkpoint_every_n=500_000
)
solve(prob, ROCK2(eigen_est=eigen_estimator), callback=callback, dt=80.0, adaptive=false)
 
# Assemble full history from all snapshot files (across all crash/restart cycles)
reload_snapshots!(saved_values, "output/", "my_simulation")
```
"""
function get_simulation_callback(; saveat, print_every_n=1000,
                                   checkpoint_dir="", checkpoint_id="latest",
                                   checkpoint_every_n=0, Float_used_to_save=Float32)
 
    use_disk = !isempty(checkpoint_dir)
    use_restart_checkpoint = use_disk && checkpoint_every_n > 0
 
    saved_values = SavedValues(Float64, Array{Float_used_to_save, 3})
 
    # --- Determine snapshot file numbering (clean stale files, continue from valid ones) ---
    snapshot_counter = Ref(0)
    if use_disk
        mkpath(checkpoint_dir)
        snapshot_counter[] = _clean_and_count_snapshots(checkpoint_dir, checkpoint_id)
        if snapshot_counter[] > 0
            println("Found $(snapshot_counter[]) valid snapshot files, continuing numbering from there.")
        end
    end
 
    # --- Snapshot saving callback (triggered at saveat times) ---
    if use_disk
        function save_func_disk(u, t, integrator)
            u_cpu = copy(adapt(CPU(), u))
            u_save = Float_used_to_save.(u_cpu)
            t_save = Float64(t)
            snapshot_counter[] += 1
            @save _snapshot_path(checkpoint_dir, checkpoint_id, snapshot_counter[]) u_save t_save
            return u_save
        end
        save_cb = SavingCallback(save_func_disk, saved_values, saveat=saveat)
    else
        function save_func_mem(u, t, integrator)
            return Float_used_to_save.(copy(adapt(CPU(), u)))
        end
        save_cb = SavingCallback(save_func_mem, saved_values, saveat=saveat)
    end
 
    # --- Print callback ---
    function print_condition(u, t, integrator)
        return integrator.stats.naccept % print_every_n == 0
    end
 
    function print_affect!(integrator)
        t = integrator.t
        step = integrator.stats.naccept
        if t > 3600 * 24 * 365
            println("Step $(step), t = $(t) s, or $(round(t / 31536000, digits=4)) years")
        elseif t > 3600 * 24
            println("Step $(step), t = $(t) s, or $(round(t / 86400, digits=2)) days")
        else
            println("Step $(step), t = $(t) s, or $(round(t / 3600, digits=2)) hours")
        end
        flush(stdout)
    end
 
    print_cb = DiscreteCallback(print_condition, print_affect!, save_positions=(false, false))
 
    # --- ADI + Advection callback (must run every step) ---
    ADI_and_ADV = DiscreteCallback((u, t, integrator) -> true, ADI_and_ADV_callback!,
                                   save_positions=(false, false))
 
    # --- Restart checkpoint callback (optional, every N steps) ---
    # Ordering in CallbackSet matters:
    #   1. ADI_and_ADV  -- state is fully consistent after this
    #   2. checkpoint    -- saves the consistent state to disk
    #   3. print         -- progress output
    #   4. save          -- snapshot at saveat times
    if use_restart_checkpoint
        cp_path = _checkpoint_path(checkpoint_dir, checkpoint_id)
 
        checkpoint_condition(u, t, integrator) = integrator.stats.naccept % checkpoint_every_n == 0
 
        function checkpoint_affect!(integrator)
            u_cpu = Float64.(Array(adapt(CPU(), integrator.u)))
            t_checkpoint = Float64(integrator.t)
            @save cp_path u_cpu t_checkpoint
            println("Restart checkpoint saved at t = $(t_checkpoint) s",
                    " ($(round(t_checkpoint / 31536000, digits=4)) years)")
            flush(stdout)
        end
 
        checkpoint_cb = DiscreteCallback(checkpoint_condition, checkpoint_affect!,
                                         save_positions=(false, false))
        callback = CallbackSet(ADI_and_ADV, checkpoint_cb, print_cb, save_cb)
    else
        callback = CallbackSet(ADI_and_ADV, print_cb, save_cb)
    end
 
    return callback, saved_values
end
 
 
# =============================================================================
# Post-solve snapshot assembly
# =============================================================================
 
"""
    reload_snapshots!(saved_values, checkpoint_dir, checkpoint_id; Float_used_to_save=Float32)
 
Load all snapshot files from disk into `saved_values`, sorted by time.
 
This replaces the contents of `saved_values.t` and `saved_values.saveval` with the full
history from all snapshot files matching the given `checkpoint_id`. Call this after `solve`
to assemble the complete history across all crash/restart cycles.
 
The ODE solver clears `saved_values` at the start of each `solve()` call, so without
calling this function, `saved_values` only contains snapshots from the most recent run.
 
# Arguments
- `saved_values`: The `SavedValues` object returned by [`get_simulation_callback`](@ref)
- `checkpoint_dir`: Directory containing snapshot files
- `checkpoint_id`: Unique identifier matching the one used in `get_simulation_callback`
- `Float_used_to_save=Float32`: Floating point type matching the one used in `get_simulation_callback`
 
# Example
```julia
callback, saved_values = get_simulation_callback(
    saveat=saveat, checkpoint_dir="output/", checkpoint_id="my_sim", checkpoint_every_n=500_000)
solve(prob, ROCK2(...), callback=callback, ...)
 
# After solve, saved_values only has this run's data.
# Reload to get the full history:
reload_snapshots!(saved_values, "output/", "my_sim")
 
# Now saved_values.t and saved_values.saveval contain ALL snapshots, sorted by time.
T_final = saved_values.saveval[end]
```
"""
function reload_snapshots!(saved_values, checkpoint_dir, checkpoint_id; Float_used_to_save=Float32)
    times, arrays, count = _load_existing_snapshots(checkpoint_dir, checkpoint_id, Float_used_to_save)
 
    # Sort by time
    perm = sortperm(times)
    times = times[perm]
    arrays = arrays[perm]
 
    # Replace contents of saved_values
    resize!(saved_values.t, count)
    resize!(saved_values.saveval, count)
    for i in 1:count
        saved_values.t[i] = times[i]
        saved_values.saveval[i] = arrays[i]
    end
 
    println("Loaded $(count) snapshots from disk into saved_values",
            count > 0 ? " (t = $(round(times[1] / 31536000, digits=4)) to $(round(times[end] / 31536000, digits=4)) years)" : "")
 
    return saved_values
end
 
 
# =============================================================================
# Restart helper
# =============================================================================
 
"""
    prepare_restart(T0, tspan, saveat; checkpoint_dir, checkpoint_id="latest", backend=CPU())
 
Check for an existing checkpoint and prepare the simulation for a fresh start or a restart.
 
If a checkpoint file exists at `checkpoint_dir/checkpoint_{checkpoint_id}.jld2`, the saved
temperature field and time are loaded. The initial condition is replaced, `tspan` is adjusted
to start from the checkpoint time, and `saveat` is filtered to only include future save times.
 
If no checkpoint file exists, all arguments are returned unchanged (identity operation).
 
This function pairs with [`get_simulation_callback`](@ref): use the same `checkpoint_dir`
and `checkpoint_id` for both so that checkpoints, snapshots, and restart logic all align.
 
# Arguments
- `T0`: Fresh initial condition (e.g., from `initial_condition_thermal_gradient`)
- `tspan`: Full time span `(t_start, t_end)`
- `saveat`: Save times (any iterable of times)
- `checkpoint_dir`: Directory where checkpoint files are stored
- `checkpoint_id="latest"`: Unique identifier matching the one used in `get_simulation_callback`
- `backend=CPU()`: Computation backend to adapt the loaded array to
 
# Returns
`(T0, tspan, saveat)` -- either unchanged (no checkpoint) or updated for restart.
 
# Example
```julia
T0_fresh = initial_condition_thermal_gradient(backend, Float64, gridx, gridy, gridz;
    T_surface=10.0, gradient=0.035)
tspan_full = (0.0, 3600.0 * 24 * 365 * 20)
saveat_full = range(tspan_full..., 21)
 
T0, tspan, saveat = prepare_restart(
    T0_fresh, tspan_full, saveat_full;
    checkpoint_dir="output/",
    checkpoint_id="my_sim",
    backend=CUDABackend()
)
 
prob = ODEProblem(rhs_diffusion_z!, T0, tspan, cache)
callback, saved_values = get_simulation_callback(saveat=saveat,
    checkpoint_dir="output/", checkpoint_id="my_sim", checkpoint_every_n=500_000)
solve(prob, ROCK2(...), callback=callback, ...)
 
# Assemble full history from disk
reload_snapshots!(saved_values, "output/", "my_sim")
```
"""
function prepare_restart(T0, tspan, saveat; checkpoint_dir, checkpoint_id="latest", backend)
    cp_path = _checkpoint_path(checkpoint_dir, checkpoint_id)
 
    if !isfile(cp_path)
        println("No checkpoint found at $(cp_path), starting fresh.")
        return T0, tspan, saveat
    end
 
    @load cp_path u_cpu t_checkpoint
 
    println("Loaded restart checkpoint from $(cp_path)")
    println("  Restarting at t = $(t_checkpoint) s ($(round(t_checkpoint / 31536000, digits=4)) years)")
 
    T0_restart = adapt(backend, eltype(T0).(u_cpu))
    tspan_restart = (t_checkpoint, tspan[2])
    saveat_restart = collect(filter(t -> t > t_checkpoint, saveat))
 
    println("  tspan: $(tspan_restart)")
    println("  saveat points remaining: $(length(saveat_restart))")
 
    return T0_restart, tspan_restart, saveat_restart
end
 
