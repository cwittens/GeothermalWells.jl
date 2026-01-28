
"""
    create_advection_index_lists(backend, gridx, gridy, gridz, boreholes)

Create index lists for advection in inner and outer pipes.

Returns tuple with index lists for efficient advection kernel dispatch.
"""
function create_advection_index_lists(backend, gridx, gridy, gridz, boreholes)

    # calculate index map
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

    u_tmp = zeros(backend, eltype(gridx), countxy_inner + countxy_outer, countz)

    return (Idx_list_Inner, Idx_list_Outer, Idx_list, countxy_inner, countxy_outer, countz, u_tmp)
end


"""
    create_cache(; backend, gridx, gridy, gridz, materials, boreholes, inlet_model)

Create simulation cache with precomputed data and temporary arrays.

Returns named tuple containing grids, materials, index lists, outlet temperature arrays,
eigenvalue estimates, and precomputed `Val` types for kernel dispatch.
"""
function create_cache(; backend, gridx, gridy, gridz, materials, boreholes, inlet_model)

    Nx, Ny, Nz = length(gridx), length(gridy), length(gridz)
    N_bh = length(boreholes)

    Idx_list_Inner, Idx_list_Outer, Idx_list, countxy_inner, countxy_outer, countz, u_tmp = create_advection_index_lists(backend, gridx, gridy, gridz, boreholes)

    T_outlet = zeros(backend, eltype(gridx), N_bh)
    T_outlet_counter = zeros(backend, Int, N_bh)

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

    # TODO: check if h is in gridz!

    gridx = adapt(backend, gridx)
    gridy = adapt(backend, gridy)
    gridz = adapt(backend, gridz)

    cache = (;
        backend,
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
        T_outlet_counter,
        u_tmp,
        Idx_list_Inner,
        Idx_list_Outer,
        Idx_list,
        countxy_inner,
        countxy_outer,
        countz,
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



# callbacks
# TODO: save T_outlet in the future


"""
    save_and_print_callback(saveat; print_every_n=1000, write_to_jld=false, data_folder_dir="", Float_used_to_save=Float32)

Create callbacks for printing progress and saving solution snapshots.

Returns `(save_cb, print_cb, saved_values)` for use with ODE solver.
Optionally writes checkpoints to JLD2 files if `write_to_jld=true`.
"""
function save_and_print_callback(saveat; print_every_n=1000, write_to_jld=false, data_folder_dir="", Float_used_to_save=Float32)

    function print_condition_open(u, t, integrator, print_every_n)
        return integrator.stats.naccept % print_every_n == 0
    end

    print_condition(u, t, integrator) = print_condition_open(u, t, integrator, print_every_n)

    function print_affect!(integrator)
        if integrator.t > 3600 * 24 * 365
            println("Step $(integrator.stats.naccept), t = $(integrator.t)s, or $(round((integrator.t / 31536000), digits=4)) years")
        elseif integrator.t > 3600 * 24
            println("Step $(integrator.stats.naccept), t = $(integrator.t)s, or $(round((integrator.t / 86400), digits=2)) days")
        else integrator.t > 3600
            println("Step $(integrator.stats.naccept), t = $(integrator.t)s, or $(round((integrator.t / 3600), digits=2)) hours")
        end
        flush(stdout)
    end

    # to have some process for long simulations
    print_cb = DiscreteCallback(print_condition, print_affect!, save_positions=(false, false))

    saved_values = SavedValues(Float64, Array{Float_used_to_save,3})

    if write_to_jld

        file_counter = Ref(0)
        if !isdir(data_folder_dir)
            mkdir(data_folder_dir)
        end

        function save_julia_array_and_write_to_JLD_open(u, t, integrator, data_folder_dir)
            u_cpu = copy(adapt(CPU(), u))
            file_counter[] += 1
            @save joinpath(data_folder_dir, "checkpoint_$(file_counter[]).jld2") u_cpu
            return u_cpu
        end

        # closer of the other function to make it work with Callback Interface
        save_julia_array_and_write_to_JLD(u, t, integrator) = save_julia_array_and_write_to_JLD_open(u, t, integrator, data_folder_dir)

        save_cb = SavingCallback(save_julia_array_and_write_to_JLD, saved_values, saveat=saveat)
    else
        function save_julia_array(u, t, integrator)
            return copy(adapt(CPU(), u))
        end

        save_cb = SavingCallback(save_julia_array, saved_values, saveat=saveat)
    end

    return save_cb, print_cb, saved_values

end


"""
    get_simulation_callback(; saveat, print_every_n=1000, write_to_jld=false, data_folder_dir="")

Create the required callback set for the simulation.

!!! warning "Required"
    This callback is essential for the simulation. It performs the ADI (Alternating Direction
    Implicit) method for horizontal diffusion and the semi-Lagrangian advection. Without this
    callback, only vertical diffusion (handled by ROCK2) is computed.

The callback combines three components:
1. **ADI + Advection**: Horizontal diffusion and fluid advection (runs every timestep)
2. **Progress printing**: Prints simulation progress every `print_every_n` steps
3. **Solution saving**: Saves temperature field at times specified by `saveat`

# Arguments
- `saveat`: Times at which to save the solution (e.g., `range(0, 3600, 10)` or `[0.0, 3600.0]`)
- `print_every_n=1000`: Print progress every N accepted timesteps
- `write_to_jld=false`: If `true`, also write checkpoints to JLD2 files
- `data_folder_dir=""`: Directory for JLD2 checkpoint files (required if `write_to_jld=true`)

# Returns
- `callback`: Combined `CallbackSet` to pass to `solve(..., callback=callback)`
- `saved_values`: `SavedValues` object containing saved solutions accessible via `saved_values.saveval`

# Example
```julia
callback, saved_values = get_simulation_callback(
    saveat=[0.0, 3600.0, 7200.0],
    print_every_n=100
)
solve(prob, ROCK2(eigen_est=eigen_estimator), callback=callback, dt=80.0, adaptive=false)

# Access saved solutions
T_final = saved_values.saveval[end]
```
"""
function get_simulation_callback(; saveat, print_every_n=1000, write_to_jld=false, data_folder_dir="")
    save_cb, print_cb, saved_values = save_and_print_callback(saveat; print_every_n=print_every_n, write_to_jld=write_to_jld, data_folder_dir=data_folder_dir)

    ADI_and_ADV = DiscreteCallback((u, t, integrator) -> true, ADI_and_ADV_callback!, save_positions=(false, false))

    callback = CallbackSet(ADI_and_ADV, print_cb, save_cb)

    return callback, saved_values
end




# TODO: semi discretize function stuff