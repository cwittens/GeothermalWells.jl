# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

"""
    rhs_diffusion_z!(dϕ, ϕ, cache, t)

Compute right-hand side for vertical (z-direction) heat diffusion.

This is the RHS function passed to OrdinaryDiffEq.jl and is integrated using an explicit 
stabilized ROCK method. Vertical diffusion is treated separately from horizontal (x,y) diffusion 
because the ADI scheme is only unconditionally stable in 2D. Since `Δz` is much larger than the 
smallest `Δx` and `Δy` (due to the fine grid resolution needed near the borehole), the explicit 
ROCK method is sufficient for the z-direction without imposing prohibitive time step restrictions.

Horizontal diffusion via the ADI method and advection are handled separately using the callback 
functionality (`ADI_and_ADV_callback!`)—a workaround to implement operator splitting within the 
OrdinaryDiffEq.jl framework.
"""
function rhs_diffusion_z!(dϕ, ϕ, cache, t)
    (; backend, gridx, gridy, gridz, Nx, Ny, Nz, Val_in_z, ValFalse, material_accessor) = cache

    # Diffusion in z-direction only. 
    # x and y directions handled using ADI, implemented as a  Diff eq callback
    diffusion_1D!(backend)(dϕ, ϕ, material_accessor, gridx, gridy, gridz, 0, Val_in_z, ValFalse, ndrange=(Nz, Ny, Nx))

    return nothing
end


@kernel inbounds = true function diffusion_1D!(dϕ, @Const(ϕ), material_accessor, @Const(gridx), @Const(gridy), @Const(gridz), dt, direction::Val{xyz}, plus_I::Val{plus_I_bool}) where {xyz,plus_I_bool}
    k, j, i = @index(Global, NTuple)

    @uniform half = eltype(ϕ)(0.5)

    rho_c = lookup_volumetric_heat_capacity(material_accessor, i, j, k, gridx, gridy, gridz)
    k_center = lookup_thermal_conductivity(material_accessor, i, j, k, gridx, gridy, gridz)

    idx_plus, idx_minus, Δ_plus, Δ_minus, k_plus, k_minus = idx_and_Δ_and_k_helper(i, j, k, gridx, gridy, gridz, material_accessor, direction)

    ϕ_kji = ϕ[k, j, i]

    # Compute 1D diffusion
    dϕ_val = @fastmath (((k_plus + k_center) * half * (ϕ[idx_plus...] - ϕ_kji) / Δ_plus
                         -
                         (k_center + k_minus) * half * (ϕ_kji - ϕ[idx_minus...]) / Δ_minus
    ) / ((Δ_plus + Δ_minus) * half)) / rho_c


    if plus_I_bool # compute (I + dt*A_i) * ϕ
        dϕ[k, j, i] = ϕ_kji + dt * dϕ_val
    else # compute A_i * ϕ
        dϕ[k, j, i] = dϕ_val
    end

    # Dirichlet-like BC: no vertical diffusion at the surface
    if xyz == :z && k == 1
        dϕ[k, j, i] = zero(eltype(ϕ))
    end
end


# give the correct indices, Δx and diffusion coefficients in the chosen direction
@inline function idx_and_Δ_and_k_helper(i, j, k, gridx, gridy, gridz, material_accessor, ::Val{:x})
    nx = length(gridx)

    # for von_Neumann BCs
    i_plus = (i == nx) ? i : i + 1
    i_minus = (i == 1) ? i : i - 1
    idx_plus = (k, j, i_plus)
    idx_minus = (k, j, i_minus)

    Δ_plus = (i == nx) ? (gridx[i] - gridx[i-1]) : (gridx[i+1] - gridx[i])
    Δ_minus = (i == 1) ? Δ_plus : (gridx[i] - gridx[i-1])

    k_plus = lookup_thermal_conductivity(material_accessor, i_plus, j, k, gridx, gridy, gridz)
    k_minus = lookup_thermal_conductivity(material_accessor, i_minus, j, k, gridx, gridy, gridz)

    return idx_plus, idx_minus, Δ_plus, Δ_minus, k_plus, k_minus
end

@inline function idx_and_Δ_and_k_helper(i, j, k, gridx, gridy, gridz, material_accessor, ::Val{:y})
    ny = length(gridy)

    # for von_Neumann BCs
    j_plus = (j == ny) ? j : j + 1
    j_minus = (j == 1) ? j : j - 1
    idx_plus = (k, j_plus, i)
    idx_minus = (k, j_minus, i)

    Δ_plus = (j == ny) ? (gridy[j] - gridy[j-1]) : (gridy[j+1] - gridy[j])
    Δ_minus = (j == 1) ? Δ_plus : (gridy[j] - gridy[j-1])

    k_plus = lookup_thermal_conductivity(material_accessor, i, j_plus, k, gridx, gridy, gridz)
    k_minus = lookup_thermal_conductivity(material_accessor, i, j_minus, k, gridx, gridy, gridz)

    return idx_plus, idx_minus, Δ_plus, Δ_minus, k_plus, k_minus
end

@inline function idx_and_Δ_and_k_helper(i, j, k, gridx, gridy, gridz, material_accessor, ::Val{:z})
    nz = length(gridz)

    # for von_Neumann BCs
    k_plus_idx = (k == nz) ? k : k + 1
    k_minus_idx = (k == 1) ? k : k - 1
    idx_plus = (k_plus_idx, j, i)
    idx_minus = (k_minus_idx, j, i)

    Δ_plus = (k == nz) ? (gridz[k] - gridz[k-1]) : (gridz[k+1] - gridz[k])
    Δ_minus = (k == 1) ? Δ_plus : (gridz[k] - gridz[k-1])

    k_plus = lookup_thermal_conductivity(material_accessor, i, j, k_plus_idx, gridx, gridy, gridz)
    k_minus = lookup_thermal_conductivity(material_accessor, i, j, k_minus_idx, gridx, gridy, gridz)

    return idx_plus, idx_minus, Δ_plus, Δ_minus, k_plus, k_minus
end

@kernel inbounds = true function thomas_I_minus_A!(U, @Const(RHS), material_accessor, @Const(gridx), @Const(gridy), @Const(gridz), dt, ::Val{N}, direction::Val{xy}) where {N,xy}
    k, ij = @index(Global, NTuple)

    @uniform Float_used = eltype(RHS)
    @uniform half = Float_used(0.5)
    @uniform one = 1

    # private memory
    b = @private Float_used (N,) # Ax = b <- this b
    lower = @private Float_used (N,) # subdiagonal of A / solution vector x later
    diagonal = @private Float_used (N,) # main diagonal of A
    upper = @private Float_used (N - 1,) # superdiagonal of A



    if direction == Val(:x)
        grid = gridx
        # load RHS into b
        for l in 1:N
            b[l] = RHS[k, ij, l]
        end
    elseif direction == Val(:y)
        grid = gridy
        # load RHS into b
        for l in 1:N
            b[l] = RHS[k, l, ij]
        end
    else
        error("Invalid direction chosen")
    end


    # Build Matrix A

    # Left boundary (l=1)
    Δ_plus = grid[2] - grid[1]
    if direction == Val(:x)
        rho_c = lookup_volumetric_heat_capacity(material_accessor, 1, ij, k, gridx, gridy, gridz)
        k_center = lookup_thermal_conductivity(material_accessor, 1, ij, k, gridx, gridy, gridz)
        k_plus = lookup_thermal_conductivity(material_accessor, 2, ij, k, gridx, gridy, gridz)
    elseif direction == Val(:y)
        rho_c = lookup_volumetric_heat_capacity(material_accessor, ij, 1, k, gridx, gridy, gridz)
        k_center = lookup_thermal_conductivity(material_accessor, ij, 1, k, gridx, gridy, gridz)
        k_plus = lookup_thermal_conductivity(material_accessor, ij, 2, k, gridx, gridy, gridz)
    else
        error("Invalid direction chosen")
    end


    # -A_i + I -> minus sign
    factor_left = -dt * half * ((k_center + k_plus) / rho_c) / (Δ_plus)^2

    diagonal[1] = -factor_left + one
    upper[1] = factor_left

    for l in 2:N-1
        Δ_minus = grid[l] - grid[l-1]
        Δ_plus = grid[l+1] - grid[l]
        # -A_i + I -> minus sign
        factor = -dt * 2 / (Δ_minus + Δ_plus)

        if direction == Val(:x)
            rho_c = lookup_volumetric_heat_capacity(material_accessor, l, ij, k, gridx, gridy, gridz)
            k_minus = lookup_thermal_conductivity(material_accessor, l - 1, ij, k, gridx, gridy, gridz)
            k_center = lookup_thermal_conductivity(material_accessor, l, ij, k, gridx, gridy, gridz)
            k_plus = lookup_thermal_conductivity(material_accessor, l + 1, ij, k, gridx, gridy, gridz)
        else # direction == Val(:y)
            rho_c = lookup_volumetric_heat_capacity(material_accessor, ij, l, k, gridx, gridy, gridz)
            k_minus = lookup_thermal_conductivity(material_accessor, ij, l - 1, k, gridx, gridy, gridz)
            k_center = lookup_thermal_conductivity(material_accessor, ij, l, k, gridx, gridy, gridz)
            k_plus = lookup_thermal_conductivity(material_accessor, ij, l + 1, k, gridx, gridy, gridz)
        end


        lower[l-1] = factor * half * (k_center + k_minus) / (Δ_minus * rho_c)
        diagonal[l] = -factor * (half * (k_center + k_plus) / (Δ_plus * rho_c) + half * (k_center + k_minus) / (Δ_minus * rho_c)) + one
        upper[l] = factor * half * (k_center + k_plus) / (Δ_plus * rho_c)
    end

    # Right boundary (l=N)
    Δ_minus = grid[N] - grid[N-1]
    if direction == Val(:x)
        rho_c = lookup_volumetric_heat_capacity(material_accessor, N, ij, k, gridx, gridy, gridz)
        k_minus = lookup_thermal_conductivity(material_accessor, N - 1, ij, k, gridx, gridy, gridz)
        k_center = lookup_thermal_conductivity(material_accessor, N, ij, k, gridx, gridy, gridz)
    else # direction == Val(:y)
        rho_c = lookup_volumetric_heat_capacity(material_accessor, ij, N, k, gridx, gridy, gridz)
        k_minus = lookup_thermal_conductivity(material_accessor, ij, N - 1, k, gridx, gridy, gridz)
        k_center = lookup_thermal_conductivity(material_accessor, ij, N, k, gridx, gridy, gridz)
    end
    # -A_i + I -> minus sign
    factor_right = -dt * half * ((k_center + k_minus) / rho_c) / (Δ_minus)^2

    lower[N-1] = factor_right
    diagonal[N] = -factor_right + one

    # thomas algorithm:
    for l in 2:N
        w = lower[l-1] / diagonal[l-1]
        diagonal[l] -= w * upper[l-1]
        b[l] -= w * b[l-1]
    end

    # 'lower' is now the cache for the solution vector
    lower[N] = b[N] / diagonal[N]
    for l in (N-1):-1:1
        lower[l] = (b[l] - upper[l] * lower[l+1]) / diagonal[l]
    end


    # copy solution back to U
    if direction == Val(:x)
        for l in 1:N
            U[k, ij, l] = lower[l]
        end

    else # direction == Val(:y)
        for l in 1:N
            U[k, l, ij] = lower[l]
        end
    end
end

"""
    advection!(ϕ, dt, t, cache, boreholes)

Apply advective heat transport in the borehole pipes using a semi-Lagrangian method.

The water flow itself is not simulated—instead, the fluid is assumed to move at prescribed 
constant velocities: `v_inner` (downward in the inner pipe) and `v_outer` (upward in the 
outer annulus). Only the temperature field is advected according to these fixed velocity profiles.

A semi-Lagrangian approach is used rather than standard explicit advection schemes because 
the high fluid velocities would impose prohibitively small time steps under the CFL constraint. 
The semi-Lagrangian method traces characteristic lines backward in time to find the departure 
point, then interpolates the temperature there using linear interpolation between grid points.
At the turnaround point at the bottom of the borehole (depth `h`), where water transitions 
from the inner pipe to the outer annulus, perfect mixing of temperature is assumed.
"""
@inline function advection!(ϕ, dt, t, cache, boreholes)
    (; u_tmp, Idx_list, Idx_list_Inner, Idx_list_Outer, count_outer_per_bh, countxy_inner, countxy_outer, countz, gridx, gridy, gridz, backend, inlet_model, T_outlet, T_outlet_counter, T_turnaround_mean) = cache

    fill!(T_outlet, 0)
    fill!(T_outlet_counter, 0)
    fill!(T_turnaround_mean, 0)

    kernel_accumulate_outlet!(backend)(T_outlet, T_outlet_counter, ϕ, Idx_list_Inner, gridz, boreholes, dt, ndrange=(countxy_inner))
    T_outlet ./= T_outlet_counter


    kernel_accumulate_turnaround_mean!(backend)(T_turnaround_mean, ϕ, Idx_list_Outer, gridz, count_outer_per_bh, boreholes, ndrange=(countz, countxy_outer))

    kernel_advection!(backend)(u_tmp, ϕ, gridx, gridy, gridz, Idx_list, Idx_list_Outer, T_turnaround_mean, countxy_inner, dt, t, boreholes, inlet_model, T_outlet, ndrange=(countz, countxy_inner + countxy_outer))

    kernel_copy_advection!(backend)(ϕ, u_tmp, Idx_list, ndrange=(countz, countxy_inner + countxy_outer))

    return nothing
end

@kernel inbounds = true function kernel_accumulate_turnaround_mean!(T_turnaround_mean, @Const(ϕ), @Const(Idx_list_Outer), @Const(gridz), @Const(count_outer_per_bh), boreholes)
    k, ij_xy = @index(Global, NTuple)

    i, j, n_bh = Idx_list_Outer[ij_xy]
    h = boreholes[n_bh].h

    if gridz[k] <= h
        # FIXME this currently assumes a uniform gird in x and y direction for the mean!
        @atomic T_turnaround_mean[k, n_bh] += (ϕ[k, j, i] / count_outer_per_bh[n_bh])
    end
end


@kernel function kernel_accumulate_outlet!(T_sum, T_outlet_counter, @Const(ϕ), @Const(Idx_list_Inner),
    @Const(gridz), boreholes, dt)
    ij_xy = @index(Global)

    i, j, n_bh = Idx_list_Inner[ij_xy]
    bh = boreholes[n_bh]
    z_max = bh.v_inner * dt

    # Loop over z in outlet region
    for (k, z) in enumerate(gridz)
        if z <= z_max
            @atomic T_sum[n_bh] += ϕ[k, j, i]
            @atomic T_outlet_counter[n_bh] += 1
        end
    end
end


@kernel inbounds = true function kernel_advection!(u_tmp, @Const(ϕ), @Const(gridx), @Const(gridy), @Const(gridz), @Const(Idx_list), @Const(Idx_list_Outer), @Const(T_turnaround_mean), countxy_inner, Δt, t, boreholes, inlet_model, T_outlet)
    k, ij_xy = @index(Global, NTuple)

    i, j, n_bh = Idx_list[ij_xy]
    x, y, z = gridx[i], gridy[j], gridz[k]

    v_inner = boreholes[n_bh].v_inner
    v_outer = boreholes[n_bh].v_outer
    h = boreholes[n_bh].h


    # Below this borehole's pipe: no advection, just preserve original
    # this only comes into play if there are different borehole heights
    # Hack to have easier indexing. (See generate cache)
    if z > h
        u_tmp[ij_xy, k] = ϕ[k, j, i]

    else
        if ij_xy <= countxy_inner # r < r_inner && z <= h
            z_departure = z + v_inner * Δt
            if z_departure > h
                # time to h
                Δt1 = (h - z) / v_inner
                # remaining time
                Δt2 = Δt - Δt1
                z_departure2 = h - v_outer * Δt2

                # FIXME if z_departure2 is close to h, k_departure_right may be in the hot rock region!
                # this is currently only fixed by added h to gridz when creating the grid
                k_departure_left, k_departure_right, α = interpolation_helper(gridz, z_departure2)

                # use mean temperature at turnaround => avoids artificial heat source from accidentally taking points from the pipe wall
                # physically this assumes perfect mixing at the turnaround (which seems justifiable)
                # Inner pipe turnaround - mean temperature from outer pipe


                u_tmp[ij_xy, k] = (1 - α) * T_turnaround_mean[k_departure_left, n_bh] + α * T_turnaround_mean[k_departure_right, n_bh]


            else
                k_departure_left, k_departure_right, α = interpolation_helper(gridz, z_departure)

                u_tmp[ij_xy, k] = (1 - α) * ϕ[k_departure_left, j, i] + α * ϕ[k_departure_right, j, i]

            end

        else # r_inner + t_inner <= r < r_outer_thickness && z <= h
            z_departure = z - v_outer * Δt
            if z_departure <= 0.0
                u_tmp[ij_xy, k] = inlet_model(n_bh, T_outlet, t)
            else
                k_departure_left, k_departure_right, α = interpolation_helper(gridz, z_departure)

                u_tmp[ij_xy, k] = (1 - α) * ϕ[k_departure_left, j, i] + α * ϕ[k_departure_right, j, i]

            end

        end
    end
end

@inline function interpolation_helper(grid, departure)
    i_departure_right = gpu_searchsortedfirst(grid, departure, 1, length(grid))
    i_departure_left = i_departure_right - 1

    # Clamp to valid indices
    i_departure_left = max(i_departure_left, 1)
    i_departure_right = min(i_departure_right, length(grid))


    x_left = grid[i_departure_left]
    x_right = grid[i_departure_right]
    α = (departure - x_left) / (x_right - x_left)

    return i_departure_left, i_departure_right, α
end


@inline function gpu_searchsortedfirst(arr, x, lo, hi)
    while lo < hi
        mid = lo + (hi - lo) ÷ 2
        @inbounds if arr[mid] < x
            lo = mid + 1
        else
            hi = mid
        end
    end
    return lo
end

@kernel inbounds = true function kernel_copy_advection!(ϕ, @Const(u_tmp), @Const(IDX_LIST))
    k, ij_xy = @index(Global, NTuple)
    i, j = IDX_LIST[ij_xy]

    ϕ[k, j, i] = u_tmp[ij_xy, k]
end


"""
    ADI_and_ADV_callback!(integrator)

Callback implementing the ADI + advection operator as part of a Strang splitting scheme.

The overall time integration uses Strang splitting between two operators:
- **Operator A**: Horizontal (x,y) diffusion via ADI + semi-Lagrangian advection (this callback)
- **Operator B**: Vertical (z) diffusion via ROCK2 (the ODE right-hand side `rhs_diffusion_z!`)

Each ROCK2 step advances by `Δt`. After each step, this callback applies operator A
by calling [`ADI_and_ADV_step!`](@ref) twice, each with `Δt/2`. This produces the
merged interior of the Strang splitting:

```
A(Δt/2) B(Δt) [A(Δt/2) A(Δt/2)] B(Δt) [A(Δt/2) A(Δt/2)] B(Δt) A(Δt/2)
                \\______  ______/         \\______  ______/
                       \\/                        \\/
                    A(Δt/2) x2 per callback = effectively A(Δt)
```

The first and last half-steps of the true Strang splitting are omitted, which introduces
a one-time first-order error that should be negligible over millions of time steps.

"""
function ADI_and_ADV_callback!(integrator)

    t = integrator.t
    Δt_half = (integrator.t - integrator.tprev) / 2

    ADI_and_ADV_step!(integrator, t, Δt_half)
    ADI_and_ADV_step!(integrator, t + Δt_half, Δt_half)

    return nothing
end


"""
    ADI_and_ADV_step!(integrator, t, Δt)

Perform one ADI + advection sub-step of size `Δt` starting at time `t`.

This executes the standard Peaceman-Rachford ADI scheme for horizontal diffusion,
interleaved with semi-Lagrangian advection. The sub-step `Δt` is itself split into
two ADI half-steps with alternating implicit directions:

**First half-step (`Δt/2`):**
1. Explicit y-diffusion: `temp = (I + Δt/2 · Aᵧ) · ϕ`
2. Semi-Lagrangian advection applied to `temp`
3. Implicit x-solve: `(I - Δt/2 · Aₓ) · ϕ = temp`

**Second half-step (`Δt/2`):**
1. Explicit x-diffusion: `temp = (I + Δt/2 · Aₓ) · ϕ`
2. Semi-Lagrangian advection applied to `temp`
3. Implicit y-solve: `(I - Δt/2 · Aᵧ) · ϕ = temp`

The advection is placed after the explicit diffusion step and before the implicit Thomas solve. 
Empirically, this is the only placement that avoids numerical instabilities. The suspected reason 
is that advection introduces large thermal gradients, and applying it before an explicit Euler-like 
step causes instability, whereas the subsequent implicit solve can handle these gradients stably. 
However, the exact theoretical justification remains uncertain.

Vertical (z) diffusion is handled separately by the main ODE right-hand side (`rhs_diffusion_z!`) 
because the ADI scheme is only unconditionally stable in 2D. However, since `Δz` is much larger 
than the smallest `Δx` and `Δy` (due to the fine grid resolution needed near the borehole), an 
explicit stabilized method (ROCK2) is sufficient for the z-direction without imposing 
prohibitive time step restrictions.

The implicit solves use the Thomas algorithm for the resulting tridiagonal systems.

# Arguments
- `integrator`: OrdinaryDiffEq integrator (provides `u`, `uprev` as working arrays, and `p` as the cache)
- `t`: Current physical time [s] at the start of this sub-step
- `Δt`: Sub-step size [s]
"""
function ADI_and_ADV_step!(integrator, t, Δt)

    ϕ = integrator.u
    temp = integrator.uprev

    (; backend, material_accessor, gridx, gridy, gridz, Nx, Ny, Nz, boreholes,
        Val_in_x, Val_in_y,
        ValTrue,
        ValNx, ValNy) = integrator.p

    ## ADI dt/2 with advection dt/2 step
    # Y direction explicit / (I + 0.5dt*A_y) * ϕ
    diffusion_1D!(backend)(temp, ϕ, material_accessor, gridx, gridy, gridz, Δt / 2, Val_in_y, ValTrue, ndrange=(Nz, Ny, Nx))

    # Advection for dt/2
    advection!(temp, Δt / 2, t, integrator.p, boreholes)

    # X direction implicit (I - 0.5dt *  A_x) \ temp
    thomas_I_minus_A!(backend)(ϕ, temp, material_accessor, gridx, gridy, gridz, Δt / 2, ValNx, Val_in_x, ndrange=(Nz, Ny))


    ## ADI dt/2 with advection dt/2 step
    # X direction explicit / (I + 0.5dt*A_x) * ϕ
    diffusion_1D!(backend)(temp, ϕ, material_accessor, gridx, gridy, gridz, Δt / 2, Val_in_x, ValTrue, ndrange=(Nz, Ny, Nx))

    # Advection for dt/2
    advection!(temp, Δt / 2, t + Δt / 2, integrator.p, boreholes)

    # Y direction implicit (I - 0.5dt *  A_y) \ temp
    thomas_I_minus_A!(backend)(ϕ, temp, material_accessor, gridx, gridy, gridz, Δt / 2, ValNy, Val_in_y, ndrange=(Nz, Nx))

    return nothing
end

"""
    eigen_estimator_pre_calculation(gridz, materials) -> λ_max

Pre-calculate eigenvalue estimate [s⁻¹] for ROCK2/ROCK4 time stepping.

Estimates maximum eigenvalue based on finest grid spacing and maximum diffusivity.
Returns `λ_max ≈ 8 * d_max / Δ_min²` where conservative factor of 8 accounts for 3D diffusion operator.
"""
function eigen_estimator_pre_calculation(gridz, materials)
    # Δx_min = minimum(diff(integrator.p.gridx))
    # Δy_min = minimum(diff(integrator.p.gridy))
    Δz_min = minimum(diff(gridz))
    Δ_min = Δz_min


    # Find maximum diffusivity across all materials
    d_max = eigen_estimator_get_dmax(materials)

    # TODO 6 comes from 3 dimensions, but here we only diffuse in z direction each step
    # maybe this can be optimized
    # Conservative estimate (factor of 8 instead of 6 for safety)
    return 8 * d_max / Δ_min^2
end

"""
    eigen_estimator(integrator)

Eigenvalue estimator function for adaptive time stepping in ROCK methods.
"""
@inline eigen_estimator(integrator) = integrator.eigen_est = integrator.p.eigen_estimate

end # end of @muladd block