using TestItems

@testitem "Directory helpers" begin
    using GeothermalWells
    @test isdir(data_dir())
    @test isdir(pkg_dir())
    @test isdir(examples_dir())
end

@testitem "data_brown_single_well_b" begin
    using GeothermalWells

    (z_beier, T_beier), (z_brown, T_brown) = data_brown_single_well_b()

    @test length(z_beier) == 118
    @test length(T_beier) == 118
    # First inlet row: T=2.3417..., z=2.3917...
    @test T_beier[1] ≈ 2.3417776680390263
    @test z_beier[1] ≈ 2.3917984132453256
    # Last outlet row (reversed, so last entry is first outlet row): T=4.7348..., z=4.0352...
    @test T_beier[end] ≈ 4.7348168668139685
    @test z_beier[end] ≈ 4.035276492085025

    @test length(z_brown) == 122
    @test length(T_brown) == 122
    # First inlet row: T=2.2404..., z=12.2088...
    @test T_brown[1] ≈ 2.2404302862412977
    @test z_brown[1] ≈ 12.208840804181094
    # Last outlet row (reversed): T=4.6679..., z=24.2500...
    @test T_brown[end] ≈ 4.667941383586896
    @test z_brown[end] ≈ 24.250056861813256
end

@testitem "data_brown_single_well_c" begin
    using GeothermalWells

    # i=1: single_300m.csv — 88 rows
    x1, T1 = data_brown_single_well_c(1)
    @test length(x1) == 88
    @test x1[1] ≈ -198.19371372535554
    @test T1[1] ≈ 19.04094909480255
    @test x1[end] ≈ 194.58592100915365
    @test T1[end] ≈ 19.038262835669904

    # i=2: single_600m.csv — 110 rows
    x2, T2 = data_brown_single_well_c(2)
    @test length(x2) == 110
    @test x2[1] ≈ -198.19371372535554
    @test T2[1] ≈ 29.06606817783625

    # i=3: single_920m.csv — 232 rows
    x3, T3 = data_brown_single_well_c(3)
    @test length(x3) == 232
    @test x3[1] ≈ -190.6869641014245
    @test T3[1] ≈ 39.67813488135304
end

@testitem "data_hu" begin
    using GeothermalWells

    # 1 year — 78 rows (two columns: depth, temperature)
    d1, t1 = data_hu(1)
    @test length(d1) == 78
    @test d1[1] ≈ 29.31877997926381
    @test t1[1] ≈ 31.745860927152314
    @test d1[end] ≈ 3468.697148195833
    @test t1[end] ≈ 36.031537978774224

    # 5 year — 91 rows
    d5, t5 = data_hu(5)
    @test length(d5) == 91

    # 10 year — 85 rows
    d10, t10 = data_hu(10)
    @test length(d10) == 85

    # 25 year — 109 rows
    d25, t25 = data_hu(25)
    @test length(d25) == 109
end

@testitem "data_li" begin
    using GeothermalWells

    # i=1: Li_500m.csv — 76 rows
    x1, T1 = data_li(1)
    @test length(x1) == 76
    @test x1[1] ≈ 0.03681585
    @test T1[1] ≈ 23.4886608
    @test x1[end] ≈ 82.40104548
    @test T1[end] ≈ 29.56358299

    # i=2: Li_1000m.csv — 90 rows
    x2, T2 = data_li(2)
    @test length(x2) == 90
    @test x2[1] ≈ 0.03681585
    @test T2[1] ≈ 23.4886608

    # i=3: Li_1500m.csv — 102 rows
    x3, T3 = data_li(3)
    @test length(x3) == 102
    @test x3[1] ≈ 0.073187778
    @test T3[1] ≈ 25.39168604

    # i=4: Li_2000m.csv — 100 rows
    x4, T4 = data_li(4)
    @test length(x4) == 100
    @test x4[1] ≈ 0.057599809
    @test T4[1] ≈ 26.94890916
end

@testitem "plot_grid smoke test" begin
    using GeothermalWells
    using Plots

    gridx = collect(-5.0:1.0:5.0)
    gridy = collect(-5.0:1.0:5.0)
    bh = Borehole{Float64}(0.0, 0.0, 2000.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)
    p = plot_grid(gridx, gridy; boreholes=(bh,))
    @test p isa Plots.Plot

    # Multiple boreholes — covers the else branch (empty labels for 2nd borehole)
    bh2 = Borehole{Float64}(3.0, 0.0, 2000.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)
    p2 = plot_grid(gridx, gridy; boreholes=(bh, bh2))
    @test p2 isa Plots.Plot

    # Zoomed-in view so annotations are triggered (x_range < 50 * r_backfill)
    gridx_zoom = collect(-0.5:0.01:0.5)
    gridy_zoom = collect(-0.5:0.01:0.5)
    p3 = plot_grid(gridx_zoom, gridy_zoom; boreholes=(bh,), annotate=true)
    @test p3 isa Plots.Plot

    # subplot parameter branch
    plot(layout=(1, 2))
    p4 = plot_grid(gridx, gridy; boreholes=(bh,), subplot=1)
    @test p4 isa Plots.Plot
end

@testitem "extract_x_profile" begin
    using GeothermalWells

    gridx = collect(-10.0:1.0:10.0)   # 21 points
    gridy = collect(-10.0:1.0:10.0)   # 21 points
    gridz = collect(0.0:100.0:1000.0)  # 11 points

    # Uniform temperature = 25.0
    T = fill(25.0, length(gridz), length(gridy), length(gridx))

    # Extract half-profile at z=500, centered at (0,0)
    r_vals, T_profile = GeothermalWells.extract_x_profile(T, gridx, gridy, gridz, 500.0, 0.0, 0.0)
    @test all(T_profile .≈ 25.0)
    @test all(r_vals .> 0)
    @test length(r_vals) == 10  # x > 0: indices for 1..10

    # Full profile
    x_full, T_full = GeothermalWells.extract_x_profile(T, gridx, gridy, gridz, 500.0, 0.0, 0.0, true)
    @test length(x_full) == length(gridx)
    @test all(T_full .≈ 25.0)

    # Test interpolation between layers
    # Set layer at z=500 (index 6) to 30.0 and z=600 (index 7) to 40.0
    T[6, :, :] .= 30.0
    T[7, :, :] .= 40.0
    _, T_interp = GeothermalWells.extract_x_profile(T, gridx, gridy, gridz, 550.0, 0.0, 0.0, true)
    @test all(T_interp .≈ 35.0)

    # z_target below grid minimum — warns and uses first layer
    T_fresh = fill(42.0, length(gridz), length(gridy), length(gridx))
    _, T_below = @test_warn "below the grid minimum" GeothermalWells.extract_x_profile(
        T_fresh, gridx, gridy, gridz, -10.0, 0.0, 0.0, true)
    @test all(T_below .≈ 42.0)

    # z_target above grid maximum — warns and uses last layer
    _, T_above = @test_warn "above the grid maximum" GeothermalWells.extract_x_profile(
        T_fresh, gridx, gridy, gridz, 2000.0, 0.0, 0.0, true)
    @test all(T_above .≈ 42.0)
end

@testitem "get_temperatures_along_z_single_well" begin
    using GeothermalWells
    using Statistics

    # Mock cache with minimal fields
    bh = Borehole{Float64}(0.0, 0.0, 500.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)
    gridz = collect(0.0:100.0:1000.0)
    # Inner/Outer indices: (i, j, bh_idx) tuples
    idx_inner = [(1, 1, 1)]
    idx_outer = [(2, 2, 1)]
    cache = (boreholes=(bh,), gridz=gridz, Idx_list_Inner=idx_inner, Idx_list_Outer=idx_outer)

    # Temperature field: T[k, j, i] with linear gradient
    T = zeros(length(gridz), 3, 3)
    for k in 1:length(gridz)
        T[k, :, :] .= 10.0 + 0.03 * gridz[k]
    end

    T_inner, T_outer, gz = GeothermalWells.get_temperatures_along_z_single_well(T, cache)
    # gridz < 500: indices 1..5 (0,100,200,300,400)
    @test length(T_inner) == 5
    @test length(T_outer) == 5
    @test length(gz) == 5
    @test T_inner[1] ≈ 10.0
    @test T_outer[1] ≈ 10.0
    @test T_inner[end] ≈ 10.0 + 0.03 * 400.0  # 22.0
end
