"""
    examples_dir()

Return path to package examples directory.
"""
examples_dir() = pkgdir(GeothermalWells, "examples")::String

"""
    data_dir()

Return path to package data directory.
"""
data_dir() = pkgdir(GeothermalWells, "data")::String

"""
    pkg_dir()

Return package root directory.
"""
pkg_dir() = pkgdir(GeothermalWells)::String

"""
    data_brown_single_well_b() -> (z_beier, T_beier), (z_brown, T_brown)

Load Brown et al. single well validation data (Figure 4b).

From: Investigating scalability of deep borehole heat exchangers (doi:10.1016/j.renene.2021.01.036)
Returns depth [m] and temperature [°C] profiles for Beier and Brown models.
"""
function data_brown_single_well_b()
    # data scrapped from paper:
    # Investigating scalability of deep borehole heat exchangers: Numerical modelling of arrays with varied modes of operation
    # https://doi.org/10.1016/j.renene.2021.01.036
    # Figure 4
    path_li = joinpath(data_dir(), "Brown_et_al")
    csv_names = ["single_beier_inlet.csv", "single_beier_outlet.csv",
        "single_brown_inlet.csv", "single_brown_outlet.csv"]

    # Read all CSV files
    data = [readdlm(joinpath(path_li, name), ','; header=false) for name in csv_names]

    # Combine inlet/outlet pairs
    z_beier = vcat(data[1][:, 2], reverse(data[2][:, 2]))
    T_beier = vcat(data[1][:, 1], reverse(data[2][:, 1]))

    z_brown = vcat(data[3][:, 2], reverse(data[4][:, 2]))
    T_brown = vcat(data[3][:, 1], reverse(data[4][:, 1]))

    return (z_beier, T_beier), (z_brown, T_brown)
end

"""
    data_brown_single_well_c(i) -> x_data, T_data

Load Brown et al. single well validation data (Figure 4c).

`i=1,2,3` for 300m, 600m, 920m depths respectively.
From: doi:10.1016/j.renene.2021.01.036
"""
function data_brown_single_well_c(i)
    # data scrapped from paper:
    # Investigating scalability of deep borehole heat exchangers: Numerical modelling of arrays with varied modes of operation
    # https://doi.org/10.1016/j.renene.2021.01.036
    # Figure 4
    path_li = joinpath(data_dir(), "Brown_et_al")
    csv_names = ["single_300m.csv", "single_600m.csv", "single_920m.csv"]


    data_path = joinpath(path_li, csv_names[i])
    data = readdlm(data_path, ','; header=false)

    x_data = data[:, 1]
    T_data = data[:, 2]

    return x_data, T_data
end

"""
    data_hu(years::Int) -> depth, temperature

Load Hu et al. validation data for specified simulation years (Figure 7).

From: Numerical modeling of a coaxial borehole heat exchanger (doi:10.1016/j.renene.2019.09.141)
Returns depth [m] and temperature [°C].
"""
function data_hu(years::Int)
    # data scrapped from paper:
    # Numerical modeling of a coaxial borehole heat exchanger to exploit geothermal energy from abandoned petroleum wells in Hinton, Alberta
    # https://doi.org/10.1016/j.renene.2019.09.141
    # Figure 7
    path_data = joinpath(data_dir(), "Hu_et_al", "Hu_$(years)year.csv")
    all_data = readdlm(path_data, ','; header=false)
    depth = all_data[:, 1]
    temperature = all_data[:, 2]
    return depth, temperature
end

"""
    data_li(i) -> x_data, T_data

Load Li et al. validation data (Figure 5).

`i=1,2,3,4` for 500m, 1000m, 1500m, 2000m depths respectively.
From: Heat extraction model and characteristics of coaxial deep borehole heat exchanger
"""
function data_li(i)
    # data scrapped from paper:
    # Heat extraction model and characteristics of coaxial deep borehole heat exchanger
    # https://doi.org/10.1016/j.renene.2021.01.036
    # Figure 5
    path_li = joinpath(data_dir(), "Li_et_al")
    csv_names = ["Li_500m.csv", "Li_1000m.csv", "Li_1500m.csv", "Li_2000m.csv"]


    data_path = joinpath(path_li, csv_names[i])
    data = readdlm(data_path, ','; header=false)

    x_data = data[:, 1]
    T_data = data[:, 2]

    return x_data, T_data
end

# TODO: some plotting functions


function plot_grid(gridx, gridy;
    boreholes=nothing,
    xlims=nothing, ylims=nothing,
    size=(800, 800),
    dpi=300,
    annotate=true,
    legend=:bottomleft,
    subplot=nothing)

    # Determine plot limits if not provided
    if isnothing(xlims)
        dx = maximum(diff(gridx))
        xlims = (minimum(gridx), maximum(gridx))
    end
    if isnothing(ylims)
        dy = maximum(diff(gridy))
        ylims = (minimum(gridy), maximum(gridy))
    end

    if isnothing(subplot)
        p = plot(
            aspect_ratio=:equal,
            xlims=xlims,
            ylims=ylims,
            xlabel="x [m]",
            ylabel="y [m]",
            legend=legend,
            size=size,
            dpi=dpi,
            grid=false,
        )
        sp = 1  # default subplot index
    else
        p = current()  # get current plot
        sp = subplot
        p = plot!(p, subplot=sp, legend=legend, xlims=xlims, ylims=ylims)  # switch to specified subplot
    end

    # Plot vertical grid lines
    for x in gridx
        vline!([x], color=:lightgray, alpha=0.9, label="", subplot=sp)
    end

    # Plot horizontal grid lines
    for y in gridy
        hline!([y], color=:lightgray, alpha=0.9, label="", subplot=sp)
    end

    if !isnothing(boreholes)
        # Angular parameter for circles
        PHI = 0:0.001:2π
        for (i, borehole) in enumerate(boreholes)
            
            if i == 1
                label_inner_tube = "Inner tube"
                label_inner_pipe_wall = "Inner pipe wall"
                label_annulus_boundary = "Annulus boundary"
                label_outer_pipe_wall = "Outer pipe wall"
                label_backfill_boundary = "Backfill boundary"

            else
                label_inner_tube = ""
                label_inner_pipe_wall = ""
                label_annulus_boundary = ""
                label_outer_pipe_wall = ""
                label_backfill_boundary = ""
            end

            # Inner pipe inner radius
            x_inner = borehole.xc .+ cos.(PHI) .* borehole.r_inner
            y_inner = borehole.yc .+ sin.(PHI) .* borehole.r_inner
            plot!(x_inner, y_inner, color=1, linewidth=2, label=label_inner_tube, subplot=sp)

            # Inner pipe outer radius
            r_inner_outer = borehole.r_inner + borehole.t_inner
            x_inner_outer = borehole.xc .+ cos.(PHI) .* r_inner_outer
            y_inner_outer = borehole.yc .+ sin.(PHI) .* r_inner_outer
            plot!(x_inner_outer, y_inner_outer, color=1, linewidth=2,
                linestyle=:dash, label=label_inner_pipe_wall, subplot=sp)

            # Outer pipe inner radius (annulus outer boundary)
            x_outer = borehole.xc .+ cos.(PHI) .* borehole.r_outer
            y_outer = borehole.yc .+ sin.(PHI) .* borehole.r_outer
            plot!(x_outer, y_outer, color=2, linewidth=2, label=label_annulus_boundary, subplot=sp)

            # Outer pipe outer radius
            r_outer_outer = borehole.r_outer + borehole.t_outer
            x_outer_outer = borehole.xc .+ cos.(PHI) .* r_outer_outer
            y_outer_outer = borehole.yc .+ sin.(PHI) .* r_outer_outer
            plot!(x_outer_outer, y_outer_outer, color=2, linewidth=2,
                linestyle=:dash, label=label_outer_pipe_wall, subplot=sp)

            # Backfill radius (if different from outer pipe outer radius)
            if borehole.r_backfill > r_outer_outer
                x_backfill = borehole.xc .+ cos.(PHI) .* borehole.r_backfill
                y_backfill = borehole.yc .+ sin.(PHI) .* borehole.r_backfill
                plot!(x_backfill, y_backfill, color=3, linewidth=2,
                    linestyle=:dot, label=label_backfill_boundary, subplot=sp)
            end

            # Add text annotations if requested and if they fit in the plot
            if annotate && xlims[1] < borehole.xc < xlims[2] && ylims[1] < borehole.yc < ylims[2]
                # Check if there's enough space for annotations
                x_range = xlims[2] - xlims[1]
                y_range = ylims[2] - ylims[1]

                # Only annotate if the borehole is reasonably sized relative to the plot
                if x_range < 50 * borehole.r_backfill && y_range < 50 * borehole.r_backfill
                    annotate!(borehole.xc, borehole.yc + borehole.r_inner / 2,
                        text("Working fluid\n(inner tube)", :center, 8))
                    annotate!(borehole.xc, borehole.yc + (r_inner_outer + borehole.r_outer) / 2,
                        text("Working fluid\n(annulus)", :center, 8))

                    # Only add grout/rock label if there's space
                    if xlims[2] - borehole.xc > 1.5 * borehole.r_backfill
                        annotate!(borehole.xc + 1.2 * borehole.r_backfill,
                            borehole.yc + 1.2 * borehole.r_backfill,
                            text("Grout/Rock", :center, 8))
                    end
                end
            end
        end
    end

    return p
end


# Function to extract 1D profile along x-axis at given depth (only x > xc)
function extract_x_profile(T, gridx, gridy, gridz, z_target, xc, yc, full_profile=false)
    # Find closest y index (at borehole center)
    y_idx = argmin(abs.(gridy .- yc))

    # Interpolate in z-direction
    if z_target <= gridz[1]
        @warn "Requested z_target $z_target is below the grid minimum $(gridz[1]). Using first layer."
        T_full = T[1, y_idx, :]
    elseif z_target >= gridz[end]
        @warn "Requested z_target $z_target is above the grid maximum $(gridz[end]). Using last layer."
        # Above grid, use last layer
        T_full = T[end, y_idx, :]
    else
        # Find bracketing indices
        z_idx_lower = findlast(gridz .<= z_target)
        z_idx_upper = findfirst(gridz .>= z_target)

        if z_idx_lower == z_idx_upper
            # Exactly on a grid point
            T_full = T[z_idx_lower, y_idx, :]
        else
            # Linear interpolation between layers
            z_lower = gridz[z_idx_lower]
            z_upper = gridz[z_idx_upper]
            weight_upper = (z_target - z_lower) / (z_upper - z_lower)
            weight_lower = 1.0 - weight_upper

            T_full = weight_lower * T[z_idx_lower, y_idx, :] + weight_upper * T[z_idx_upper, y_idx, :]
        end
    end

    # Filter for x > xc
    if full_profile
        return gridx, T_full
    else
        mask = gridx .> xc
        x_filtered = gridx[mask]
        T_profile = T_full[mask]

        # Compute radial distances (now just x - xc, not absolute value)
        r_vals = x_filtered .- xc

        return r_vals, T_profile
    end
end

"""
    get_temperatures_along_z_single_well(ϕ_t, cache)

Extract average temperatures in inner and outer pipes as a function of depth.

Returns `(T_inner, T_outer, gridz_adjusted)` where temperatures are averaged
over all grid points in each pipe section at each depth level.

Note: This function is only valid for single borehole configurations.
"""
function get_temperatures_along_z_single_well(ϕ_t, cache)
    if length(cache.boreholes) != 1
        @warn "This function is only correct for single borehole setups"
    end
    
    h = cache.boreholes[1].h
    gridz = cache.gridz
    idx_inner = cache.Idx_list_Inner
    idx_outer = cache.Idx_list_Outer
    
    # Find depth indices up to borehole depth
    Nz_till_h = sum(gridz .< h)
    
    T_inner = Float64[]
    T_outer = Float64[]
    
    for k in 1:Nz_till_h
        # Collect all temperatures at this depth level
        temp_inner = [ϕ_t[k, j, i] for (i, j, bh_idx) in idx_inner]
        temp_outer = [ϕ_t[k, j, i] for (i, j, bh_idx) in idx_outer]
        
        # Average over all grid points at this depth
        push!(T_inner, mean(temp_inner))
        push!(T_outer, mean(temp_outer))
    end
    
    return T_inner, T_outer, gridz[1:Nz_till_h]
end