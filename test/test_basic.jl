using TestItems

@testitem "Package loads" begin
    using GeothermalWells
    @test true
end

@testitem "Borehole construction" begin
    using GeothermalWells

    # Based on Hu et al. - deep borehole with insulation
    bh = Borehole{Float64}(
        0, 0,        # xc, yc
        3500,          # h - 3500m depth
        0.0381,          # r_inner
        0.01,            # t_inner
        0.0889,          # r_outer
        0.01,            # t_outer
        0.0989,          # r_backfill (r_outer + t_outer)
        10,            # ṁ = 10 kg/s
        1000           # insulation_depth
    )
    # todo add values for v_inner and 
    @test bh.h == 3500.0
    @test bh.ṁ == 10.0
    @test bh.v_inner > 0
    @test bh.v_outer > 0
    @test bh.insulation_depth == 1000.0
end



@testitem "HomogenousMaterialProperties construction" begin
    using GeothermalWells

    materials = HomogenousMaterialProperties{Float64}(
        2.88,              # k_rock
        2.17e6,            # rho_c_rock
        0.6,               # k_water
        4.18e6,            # rho_c_water
        44.5,              # k_steel
        3.73e6,            # rho_c_steel
        0.26,              # k_insulating
        1.96e6,            # rho_c_insulating
        1.0,               # k_backfill (dummy)
        1.0                # rho_c_backfill (dummy)
    )

    @test materials.k_rock == 2.88
    @test materials.k_water == 0.6
end

@testitem "StratifiedMaterialProperties" begin
    using GeothermalWells

    # 4 stratified rock layers from Li et al.
    materials = StratifiedMaterialProperties{4,Float64}(
        (1.8, 2.6, 3.5, 5.3),                                    # k_rock per layer
        (1780*1379.0, 2030*1450.0, 1510*1300.0, 2600*878.0),     # rho_c_rock per layer
        (500.0, 1000.0, 1500.0, 2000.0),                         # layer depths
        0.618, 4.166e6,    # water
        41.0, 7850*475.0,  # steel
        0.4, 1.955e6,      # insulating
        1.5, 1.76e6        # backfill
    )

    @test materials.k_rock_layers[1] == 1.8
    @test materials.k_rock_layers[4] == 5.3
    @test materials.layer_depths[2] == 1000.0

    # get_rock_thermal_conductivity_by_depth — picks correct layer
    @test GeothermalWells.get_rock_thermal_conductivity_by_depth(100.0, materials) == 1.8   # layer 1 (z ≤ 500)
    @test GeothermalWells.get_rock_thermal_conductivity_by_depth(500.0, materials) == 1.8   # boundary of layer 1
    @test GeothermalWells.get_rock_thermal_conductivity_by_depth(750.0, materials) == 2.6   # layer 2 (500 < z ≤ 1000)
    @test GeothermalWells.get_rock_thermal_conductivity_by_depth(1250.0, materials) == 3.5  # layer 3
    @test GeothermalWells.get_rock_thermal_conductivity_by_depth(1750.0, materials) == 5.3  # layer 4
    @test GeothermalWells.get_rock_thermal_conductivity_by_depth(3000.0, materials) == 5.3  # below all layers → deepest

    # get_rock_volumetric_heat_capacity_by_depth — picks correct layer
    @test GeothermalWells.get_rock_volumetric_heat_capacity_by_depth(100.0, materials) == 1780 * 1379.0
    @test GeothermalWells.get_rock_volumetric_heat_capacity_by_depth(750.0, materials) == 2030 * 1450.0
    @test GeothermalWells.get_rock_volumetric_heat_capacity_by_depth(1250.0, materials) == 1510 * 1300.0
    @test GeothermalWells.get_rock_volumetric_heat_capacity_by_depth(1750.0, materials) == 2600 * 878.0
    @test GeothermalWells.get_rock_volumetric_heat_capacity_by_depth(3000.0, materials) == 2600 * 878.0

    # eigen_estimator_get_dmax — max diffusivity across all materials
    d_max = GeothermalWells.eigen_estimator_get_dmax(materials)
    d_rock_layers = [materials.k_rock_layers[i] / materials.rho_c_rock_layers[i] for i in 1:4]
    d_water = materials.k_water / materials.rho_c_water
    d_steel = materials.k_steel / materials.rho_c_steel
    d_insulating = materials.k_insulating / materials.rho_c_insulating
    d_backfill = materials.k_backfill / materials.rho_c_backfill
    @test d_max ≈ max(maximum(d_rock_layers), d_water, d_steel, d_insulating, d_backfill)
end

@testitem "Grid creation" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    bh = Borehole{Float32}(0.0, 0.0, 2000 - 1, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)

    # Test domain computation
    domain = compute_domain([bh]; buffer_x=100, buffer_y=100, buffer_z=200+1)
    @test domain.xmin == -100
    @test domain.xmax == 100
    @test domain.zmax == 2200

    # Test uniform z-grid includes borehole depth
    gridz = create_uniform_gridz_with_borehole_depths(zmin = 0, zmax = 2200, dz = 100, boreholes = (bh,), backend = CPU())
    @test eltype(gridz[1]) == Float32
    @test bh.h ∈ gridz
    @test gridz[1] == 0
end

@testitem "Adaptive grid 1D" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    bh = Borehole{Float64}(0.0, 0.0, 1999.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)

    gridx = create_adaptive_grid_1d(
        xmin=-100, xmax=100,
        dx_fine=0.0025, growth_factor=1.3, dx_max=10.0,
        boreholes=(bh,), backend=CPU(), Float_used=Float64, direction=:x
    )

    @test gridx[1] == -100
    @test gridx[end] == 100
    @test minimum(diff(gridx)) <= 0.0025
end

@testitem "Adaptive grid 1D - multiple boreholes" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    # Two boreholes separated by 30m in x
    bh1 = Borehole{Float64}(-15.0, 0.0, 2000.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)
    bh2 = Borehole{Float64}(15.0, 0.0, 2000.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)
    boreholes = (bh1, bh2)

    gridx = create_adaptive_grid_1d(
        xmin=-100, xmax=100,
        dx_fine=0.0025, growth_factor=1.3, dx_max=10.0,
        boreholes=boreholes, backend=CPU(), Float_used=Float64, direction=:x
    )

    @test gridx[1] == -100.0
    @test gridx[end] == 100.0
    @test issorted(gridx)
    @test length(unique(gridx)) == length(gridx)  # no duplicates

    # Fine spacing near both boreholes
    near_bh1 = filter(x -> abs(x - bh1.xc) < bh1.r_backfill, gridx)
    near_bh2 = filter(x -> abs(x - bh2.xc) < bh2.r_backfill, gridx)
    @test length(near_bh1) > 0
    @test length(near_bh2) > 0

    # The gap between boreholes should have grid points
    between = filter(x -> bh1.xc + bh1.r_backfill < x < bh2.xc - bh2.r_backfill, gridx)
    @test length(between) > 0

    # Also test y-direction with the same boreholes offset in y
    bh3 = Borehole{Float64}(0.0, -15.0, 2000.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)
    bh4 = Borehole{Float64}(0.0, 15.0, 2000.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)

    gridy = create_adaptive_grid_1d(
        xmin=-100, xmax=100,
        dx_fine=0.0025, growth_factor=1.3, dx_max=10.0,
        boreholes=(bh3, bh4), backend=CPU(), Float_used=Float64, direction=:y
    )

    @test gridy[1] == -100.0
    @test gridy[end] == 100.0
    @test issorted(gridy)
end

@testitem "Uniform z-grid - multiple boreholes with different depths" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    bh1 = Borehole{Float64}(-15.0, 0.0, 1500.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)
    bh2 = Borehole{Float64}(15.0, 0.0, 2000.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)

    gridz = create_uniform_gridz_with_borehole_depths(
        zmin=0, zmax=2500, dz=100,
        boreholes=(bh1, bh2), backend=CPU()
    )

    @test gridz[1] == 0.0
    @test bh1.h ∈ gridz  # 1500.0 included
    @test bh2.h ∈ gridz  # 2000.0 included
    @test issorted(gridz)
    @test length(unique(gridz)) == length(gridz)
end

@testitem "Initial condition - thermal gradient" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    gridx = collect(0.0:10.0:100.0)
    gridy = collect(0.0:10.0:100.0)
    gridz = collect(0.0:100.0:3700.0)

    # Hu et al. style: T_surface=2.29, gradient=0.035 K/m
    T0 = initial_condition_thermal_gradient(
        CPU(), Float64, gridx, gridy, gridz;
        T_surface=2.29, gradient=0.035
    )

    @test size(T0) == (length(gridz), length(gridy), length(gridx))
    @test T0[1, 1, 1] ≈ 2.29                     # surface
    @test T0[end, 1, 1] ≈ 2.29 + 0.035 * 3700.0  # bottom (~132°C)
end

@testitem "Inlet models" begin
    using GeothermalWells

    # Constant inlet (Hu et al. style)
    const_inlet = ConstantInlet{Float64}(20.0)
    @test const_inlet(1, [0.0], 0.0) == 20.0

    # Heat exchanger inlet (Li et al. style)
    # Q = 200 kW, ṁ = 11.65 kg/s, c = 4174 J/(kg·K)
    Q = 200e3
    ṁ = 42 * 998.2 / 3600  # ~11.65 kg/s
    c = 4174.0
    ΔT = Q / (ṁ * c)
    hx_inlet = HeatExchangerInlet{Float64}(ΔT)
    
    T_outlet = [50.0]
    @test hx_inlet(1, T_outlet, 0.0) ≈ 50.0 - ΔT

    # Custom inlet
    custom_inlet = CustomInlet((_,_,t) -> 20.0 + 0.01*t)
    @test custom_inlet(1, [0.0], 100.0) ≈ 21.0
end

@testitem "Multiple boreholes - array setup" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    # 3x3 array like NxM_Array.jl
    XC = [-30.0, 0.0, 30.0]
    YC = [-30.0, 0.0, 30.0]

    boreholes = tuple(
        (Borehole{Float64}(xc, yc, 2000.0, 0.0511, 0.0114, 0.0885, 0.00833, 0.115, 11.65, 0.0)
        for xc in XC, yc in YC)...
    )

    @test length(boreholes) == 9

    domain = compute_domain(boreholes; buffer_x=100, buffer_y=100, buffer_z=200)
    @test domain.xmin == -130
    @test domain.xmax == 130
    @test domain.ymin == -130
    @test domain.ymax == 130
end