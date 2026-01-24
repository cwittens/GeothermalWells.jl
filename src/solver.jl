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
    (; backend, gridx, gridy, gridz, Nx, Ny, Nz, Val_in_z, ValFalse, boreholes, materials) = cache

    # Diffusion in z-direction only. 
    # x and y directions handled using ADI, implemented as a  Diff eq callback
    diffusion_1D!(backend)(dϕ, ϕ, gridx, gridy, gridz, boreholes, materials, 0, Val_in_z, ValFalse, ndrange=(Nz, Ny, Nx))

    return nothing
end


@kernel inbounds = true function diffusion_1D!(dϕ, @Const(ϕ), @Const(gridx), @Const(gridy), @Const(gridz), boreholes, materials, dt, direction::Val{xyz}, plus_I::Val{plus_I_bool}) where {xyz,plus_I_bool}
    k, j, i = @index(Global, NTuple)

    @uniform half = eltype(ϕ)(0.5)

    rho_c = get_volumetric_heat_capacity(gridx[i], gridy[j], gridz[k], boreholes, materials)
    k_center = get_thermal_conductivity(gridx[i], gridy[j], gridz[k], boreholes, materials)

    idx_plus, idx_minus, Δ_plus, Δ_minus, k_plus, k_minus = idx_and_Δ_and_k_helper(i, j, k, gridx, gridy, gridz, boreholes, materials, direction)

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
end


# give the correct indices, Δx and diffusion coefficients in the chosen direction
@inline function idx_and_Δ_and_k_helper(i, j, k, gridx, gridy, gridz, boreholes, materials, ::Val{:x})
    nx = length(gridx)

    # for von_Neumann BCs
    idx_plus = (i == nx) ? (k, j, i) : (k, j, i + 1)
    idx_minus = (i == 1) ? (k, j, i) : (k, j, i - 1)

    Δ_plus = (i == nx) ? (gridx[i] - gridx[i-1]) : (gridx[i+1] - gridx[i])
    Δ_minus = (i == 1) ? Δ_plus : (gridx[i] - gridx[i-1])

    k_plus = get_thermal_conductivity(gridx[idx_plus[3]], gridy[j], gridz[k], boreholes, materials)
    k_minus = get_thermal_conductivity(gridx[idx_minus[3]], gridy[j], gridz[k], boreholes, materials)

    return idx_plus, idx_minus, Δ_plus, Δ_minus, k_plus, k_minus
end

@inline function idx_and_Δ_and_k_helper(i, j, k, gridx, gridy, gridz, boreholes, materials, ::Val{:y})
    ny = length(gridy)

    # for von_Neumann BCs
    idx_plus = (j == ny) ? (k, j, i) : (k, j + 1, i)
    idx_minus = (j == 1) ? (k, j, i) : (k, j - 1, i)

    Δ_plus = (j == ny) ? (gridy[j] - gridy[j-1]) : (gridy[j+1] - gridy[j])
    Δ_minus = (j == 1) ? Δ_plus : (gridy[j] - gridy[j-1])

    k_plus = get_thermal_conductivity(gridx[i], gridy[idx_plus[2]], gridz[k], boreholes, materials)
    k_minus = get_thermal_conductivity(gridx[i], gridy[idx_minus[2]], gridz[k], boreholes, materials)


    return idx_plus, idx_minus, Δ_plus, Δ_minus, k_plus, k_minus
end

@inline function idx_and_Δ_and_k_helper(i, j, k, gridx, gridy, gridz, boreholes, materials, ::Val{:z})
    nz = length(gridz)

    # for von_Neumann BCs
    idx_plus = (k == nz) ? (k, j, i) : (k + 1, j, i)
    idx_minus = (k == 1) ? (k, j, i) : (k - 1, j, i)

    Δ_plus = (k == nz) ? (gridz[k] - gridz[k-1]) : (gridz[k+1] - gridz[k])
    Δ_minus = (k == 1) ? Δ_plus : (gridz[k] - gridz[k-1])

    k_plus = get_thermal_conductivity(gridx[i], gridy[j], gridz[idx_plus[1]], boreholes, materials)
    k_minus = get_thermal_conductivity(gridx[i], gridy[j], gridz[idx_minus[1]], boreholes, materials)

    return idx_plus, idx_minus, Δ_plus, Δ_minus, k_plus, k_minus
end

@kernel inbounds = true function thomas_I_minus_A!(U, RHS, gridx, gridy, gridz, dt, boreholes, materials, ::Val{N}, direction::Val{xy}) where {N,xy}
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
        rho_c = get_volumetric_heat_capacity(gridx[1], gridy[ij], gridz[k], boreholes, materials)
        k_center = get_thermal_conductivity(gridx[1], gridy[ij], gridz[k], boreholes, materials)
        k_plus = get_thermal_conductivity(gridx[2], gridy[ij], gridz[k], boreholes, materials)
    elseif direction == Val(:y)
        rho_c = get_volumetric_heat_capacity(gridx[ij], gridy[1], gridz[k], boreholes, materials)
        k_center = get_thermal_conductivity(gridx[ij], gridy[1], gridz[k], boreholes, materials)
        k_plus = get_thermal_conductivity(gridx[ij], gridy[2], gridz[k], boreholes, materials)
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
            rho_c = get_volumetric_heat_capacity(gridx[l], gridy[ij], gridz[k], boreholes, materials)
            k_minus = get_thermal_conductivity(gridx[l-1], gridy[ij], gridz[k], boreholes, materials)
            k_center = get_thermal_conductivity(gridx[l], gridy[ij], gridz[k], boreholes, materials)
            k_plus = get_thermal_conductivity(gridx[l+1], gridy[ij], gridz[k], boreholes, materials)
        else # direction == Val(:y)
            rho_c = get_volumetric_heat_capacity(gridx[ij], gridy[l], gridz[k], boreholes, materials)
            k_minus = get_thermal_conductivity(gridx[ij], gridy[l-1], gridz[k], boreholes, materials)
            k_center = get_thermal_conductivity(gridx[ij], gridy[l], gridz[k], boreholes, materials)
            k_plus = get_thermal_conductivity(gridx[ij], gridy[l+1], gridz[k], boreholes, materials)
        end


        lower[l-1] = factor * half * (k_center + k_minus) / (Δ_minus * rho_c)
        diagonal[l] = -factor * (half * (k_center + k_plus) / (Δ_plus * rho_c) + half * (k_center + k_minus) / (Δ_minus * rho_c)) + one
        upper[l] = factor * half * (k_center + k_plus) / (Δ_plus * rho_c)
    end

    # Right boundary (l=N)
    Δ_minus = grid[N] - grid[N-1]
    if direction == Val(:x)
        rho_c = get_volumetric_heat_capacity(gridx[N], gridy[ij], gridz[k], boreholes, materials)
        k_minus = get_thermal_conductivity(gridx[N-1], gridy[ij], gridz[k], boreholes, materials)
        k_center = get_thermal_conductivity(gridx[N], gridy[ij], gridz[k], boreholes, materials)
    else # direction == Val(:y)
        rho_c = get_volumetric_heat_capacity(gridx[ij], gridy[N], gridz[k], boreholes, materials)
        k_minus = get_thermal_conductivity(gridx[ij], gridy[N-1], gridz[k], boreholes, materials)
        k_center = get_thermal_conductivity(gridx[ij], gridy[N], gridz[k], boreholes, materials)
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
    (; u_tmp, Idx_list, Idx_list_Inner, Idx_list_Outer, countxy_inner, countxy_outer, countz, gridx, gridy, gridz, backend, inlet_model, T_outlet, T_outlet_counter) = cache

    fill!(T_outlet, 0)
    fill!(T_outlet_counter, 0)

    kernel_accumulate_outlet!(backend)(T_outlet, T_outlet_counter, ϕ, Idx_list_Inner, gridz, boreholes, dt, ndrange=(countxy_inner))
    T_outlet ./= T_outlet_counter

    kernel_advection!(backend)(u_tmp, ϕ, gridx, gridy, gridz, Idx_list, Idx_list_Outer, countxy_inner, dt, t, boreholes, inlet_model, T_outlet, ndrange=(countz, countxy_inner + countxy_outer))

    kernel_copy_advection!(backend)(ϕ, u_tmp, Idx_list, ndrange=(countz, countxy_inner + countxy_outer))

    return nothing
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


@kernel inbounds = true function kernel_advection!(u_tmp, @Const(ϕ), @Const(gridx), @Const(gridy), @Const(gridz), @Const(Idx_list), @Const(Idx_list_Outer), countxy_inner, Δt, t, boreholes, inlet_model, T_outlet)
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

                # FIXME this can be made more efficient by precomputing it. 
                # use mean temperature at turnaround => avoids artificial heat source from accidentally taking points from the pipe wall
                # physically this assumes perfect mixing at the turnaround (which seems justifiable)
                # Inner pipe turnaround - mean temperature from outer pipe

                # FIXME this currently assumes a uniform gird in x and y direction for the mean!
                mean_left = 0
                mean_right = 0
                count_outer = 0

                for (i_outer, j_outer, bh_idx_outer) in Idx_list_Outer
                    if bh_idx_outer == n_bh  # Only average over THIS borehole's outer pipe
                        mean_left += ϕ[k_departure_left, j_outer, i_outer]
                        mean_right += ϕ[k_departure_right, j_outer, i_outer]
                        count_outer += 1
                    end
                end

                mean_left /= count_outer
                mean_right /= count_outer

                u_tmp[ij_xy, k] = (1 - α) * mean_left + α * mean_right


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

# TODO: better splitting like (Δt/2 ADI+ADV) + (Δt ROCK z) + (Δt/2 ADI+ADV)??
"""
    ADI_and_ADV_callback!(integrator)

Alternating Direction Implicit (ADI) callback for horizontal (x,y) diffusion combined with advection.

This callback implements operator splitting for the horizontal diffusion using the
ADI scheme, interleaved with semi-Lagrangian advection. Each full timestep `Δt` is split into 
two half-steps with alternating implicit directions:

**First half-step (Δt/2):**
1. Explicit y-diffusion: `temp = (I + Δt/2 · Aᵧ) · ϕ`
2. Advection applied to `temp`
3. Implicit x-solve: `(I - Δt/2 · Aₓ) · ϕ = temp`

**Second half-step (Δt/2):**
1. Explicit x-diffusion: `temp = (I + Δt/2 · Aₓ) · ϕ`
2. Advection applied to `temp`
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

The implicit x/y solves use the Thomas algorithm for the resulting tridiagonal systems.
"""
function ADI_and_ADV_callback!(integrator)
    Δt = integrator.t - integrator.tprev
    ϕ = integrator.u
    temp = integrator.uprev

    (; backend, gridx, gridy, gridz, Nx, Ny, Nz, boreholes, materials,
        Val_in_x, Val_in_y,
        ValTrue,
        ValNx, ValNy) = integrator.p

    ## ADI dt/2 with advection dt/2 step
    # Y direction explicit / (I + 0.5dt*A_y) * ϕ
    diffusion_1D!(backend)(temp, ϕ, gridx, gridy, gridz, boreholes, materials, Δt / 2, Val_in_y, ValTrue, ndrange=(Nz, Ny, Nx))

    # Advection for dt/2
    advection!(temp, Δt / 2, integrator.t, integrator.p, boreholes)

    # X direction implicit (I - 0.5dt *  A_x) \ temp
    thomas_I_minus_A!(backend)(ϕ, temp, gridx, gridy, gridz, Δt / 2, boreholes, materials, ValNx, Val_in_x, ndrange=(Nz, Ny))


    ## ADI dt/2 with advection dt/2 step
    # X direction explicit / (I + 0.5dt*A_x) * ϕ
    diffusion_1D!(backend)(temp, ϕ, gridx, gridy, gridz, boreholes, materials, Δt / 2, Val_in_x, ValTrue, ndrange=(Nz, Ny, Nx))

    # Advection for dt/2
    advection!(temp, Δt / 2, integrator.t + Δt / 2, integrator.p, boreholes)

    # Y direction implicit (I - 0.5dt *  A_y) \ temp
    thomas_I_minus_A!(backend)(ϕ, temp, gridx, gridy, gridz, Δt / 2, boreholes, materials, ValNy, Val_in_y, ndrange=(Nz, Nx))

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