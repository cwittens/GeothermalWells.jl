"""
    compute_domain(boreholes; buffer_x=100, buffer_y=100, buffer_z=200)

Compute simulation domain bounds from borehole positions.

Returns `(xmin, xmax, ymin, ymax, zmin, zmax)` encompassing all boreholes with specified buffers [m].
"""
function compute_domain(boreholes; buffer_x=100, buffer_y=100, buffer_z=200)
    xs = [bh.xc for bh in boreholes]
    ys = [bh.yc for bh in boreholes]
    zmax = maximum(bh.h for bh in boreholes) + buffer_z

    xmin, xmax = minimum(xs) - buffer_x, maximum(xs) + buffer_x
    ymin, ymax = minimum(ys) - buffer_y, maximum(ys) + buffer_y
    zmin = zero(zmax)

    return (; xmin, xmax, ymin, ymax, zmin, zmax)
end



"""
    create_uniform_gridz_with_borehole_depths(zmin, zmax, dz, boreholes, backend)

Create vertical grid including all borehole depths.

Ensures each borehole depth `h` is a grid point for accurate advection interpolation.
"""
function create_uniform_gridz_with_borehole_depths(zmin, zmax, dz, boreholes, backend)
    Float_used = eltype(boreholes[1].h)
    # Create uniform grid
    gridz = collect(zmin:dz:zmax)
    
    # Add all borehole depths
    # The way the advection works, when h is not part of the grid, it can currently interpolate
    # between 
    for bh in boreholes
        push!(gridz, bh.h)
    end
    
    # Remove duplicates and sort
    unique!(sort!(gridz))

    # move to backend
    gridz = adapt(backend, Float_used.(gridz))
    
    return gridz
end


"""
    create_adaptive_grid_1d(; xmin, xmax, dx_fine, growth_factor, dx_max, boreholes, backend, Float_used, direction)

Generate a non-uniform 1D grid with refined spacing around one or multiple boreholes.

Creates a 1D grid along the specified direction with fine resolution near boreholes
and geometrically increasing spacing away from them, up to a maximum cell size. Works
for both single borehole and multi-borehole configurations.

# Arguments
- `xmin`, `xmax`: Domain boundaries [m]
- `dx_fine`: Fine grid spacing near boreholes [m]
- `growth_factor`: Geometric growth rate for coarsening (e.g., 1.3)
- `dx_max`: Maximum grid spacing [m]
- `boreholes`: Single borehole or collection of borehole objects with positions and radii
- `backend`: Computation backend (CPU/GPU) for array adaptation
- `Float_used`: Floating point type for grid values
- `direction`: Grid direction, either `:x` or `:y`

# Grid generation strategy
For single or multiple boreholes:
1. Creates fine uniform grids around each borehole (spacing ≈ `dx_fine`)
2. Extends from domain boundaries toward boreholes with geometric growth
3. For multiple boreholes: fills gaps between boreholes symmetrically with progressive coarsening
4. Ensures all grid points are unique and sorted

# Returns
Backend-adapted 1D grid array of the specified floating point type.
"""
function create_adaptive_grid_1d(; xmin, xmax, dx_fine, growth_factor, dx_max, boreholes, backend, Float_used, direction)

    if direction == :x
        perm = sortperm([bh.xc for bh in boreholes])
        centers = [bh.xc for bh in boreholes][perm]
    elseif direction == :y
        perm = sortperm([bh.yc for bh in boreholes])
        centers = [bh.yc for bh in boreholes][perm]
    end

    radii = (1.0 .* [bh.r_backfill for bh in boreholes])[perm]
    centers
    # N_boreholes many intervals around each borehole with roughly dx_fine spacing
    intervalls = [collect(range(centers[i] - radii[i], centers[i] + radii[i], Int(ceil((2radii[i]) / dx_fine)) + 1)) for i in 1:length(boreholes)]


    # fill up from left most point to xmin
    dx = dx_fine
    x = intervalls[1][1]
    while x > xmin
        dx = min(dx * growth_factor, dx_max)
        x = x - dx
        if x > xmin
            pushfirst!(intervalls[1], x)
        else
            pushfirst!(intervalls[1], xmin)
        end
    end

    # fill up from right most point to xmax
    dx = dx_fine
    x = intervalls[end][end]
    while x < xmax
        dx = min(dx * growth_factor, dx_max)
        x = x + dx
        if x < xmax
            push!(intervalls[end], x)
        else
            push!(intervalls[end], xmax)
        end
    end

    # inbetween boreholes
    for i in 1:(length(boreholes)-1)
        x_i = intervalls[i][end]
        x_i_plus_1 = intervalls[i+1][1]
        dx = dx_fine
        # while x_i < x_i_plus_1
        while x_i_plus_1 - x_i > growth_factor * dx
            dx = min(dx * growth_factor, dx_max)
            x_i = x_i + dx
            x_i_plus_1 = x_i_plus_1 - dx

            if x_i_plus_1 - x_i >  dx
                push!(intervalls[i], x_i)
                pushfirst!(intervalls[i+1], x_i_plus_1)
            else
                push!(intervalls[i], (x_i + x_i_plus_1) / 2)
            end
        end
    end

    grid = unique(sort(vcat(intervalls...)))

    grid = adapt(backend, Float_used.(grid))

    return grid
end
