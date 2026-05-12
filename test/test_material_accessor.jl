using TestItems

@testitem "_material_mode dispatch" begin
    using GeothermalWells

    # :auto picks based on N_bh
    @test GeothermalWells._material_mode(:auto, 1) === :on_the_fly
    @test GeothermalWells._material_mode(:auto, 2) === :precomputed
    @test GeothermalWells._material_mode(:auto, 9) === :precomputed

    # Explicit precomputed
    @test GeothermalWells._material_mode(true, 1) === :precomputed
    @test GeothermalWells._material_mode(true, 5) === :precomputed
    @test GeothermalWells._material_mode(:precomputed, 1) === :precomputed

    # Explicit on-the-fly (all three spellings)
    @test GeothermalWells._material_mode(false, 1) === :on_the_fly
    @test GeothermalWells._material_mode(false, 5) === :on_the_fly
    @test GeothermalWells._material_mode(:on_the_fly, 5) === :on_the_fly
    @test GeothermalWells._material_mode(:onthefly, 5) === :on_the_fly

    # Invalid values
    @test_throws ArgumentError GeothermalWells._material_mode(:bogus, 1)
    @test_throws ArgumentError GeothermalWells._material_mode("precomputed", 1)
    @test_throws ArgumentError GeothermalWells._material_mode(42, 1)
end

@testitem "PrecomputedMaterialAccessor lookup uses [k, j, i] layout" begin
    using GeothermalWells

    # Build arrays with shape (Nz, Ny, Nx) = (4, 3, 2) and unique values
    Nz, Ny, Nx = 4, 3, 2
    k_arr = reshape(collect(1.0:Nz*Ny*Nx), Nz, Ny, Nx)
    rho_c_arr = reshape(collect(101.0:100.0+Nz*Ny*Nx), Nz, Ny, Nx)

    acc = GeothermalWells.PrecomputedMaterialAccessor(k_arr, rho_c_arr)
    @test acc isa GeothermalWells.AbstractMaterialAccessor

    # Grids are ignored by the precomputed accessor but the API requires them
    gridx = collect(0.0:Float64(Nx - 1))
    gridy = collect(0.0:Float64(Ny - 1))
    gridz = collect(0.0:Float64(Nz - 1))

    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        @test GeothermalWells.lookup_thermal_conductivity(acc, i, j, k, gridx, gridy, gridz) == k_arr[k, j, i]
        @test GeothermalWells.lookup_volumetric_heat_capacity(acc, i, j, k, gridx, gridy, gridz) == rho_c_arr[k, j, i]
    end
end

@testitem "OnTheFlyMaterialAccessor lookup matches get_thermal_conductivity / get_volumetric_heat_capacity" begin
    using GeothermalWells

    materials = HomogenousMaterialProperties{Float64}(
        2.88, 2.17e6, 0.6, 4.18e6, 44.5, 3.73e6, 0.26, 1.96e6, 1.0, 1.0
    )
    bh = Borehole{Float64}(0.0, 0.0, 100.0, 0.0381, 0.01, 0.0889, 0.01, 0.0989, 10.0, 50.0)
    boreholes = (bh,)

    acc = GeothermalWells.OnTheFlyMaterialAccessor(boreholes, materials)
    @test acc isa GeothermalWells.AbstractMaterialAccessor

    # Sample points covering: inner pipe, pipe walls, outer annulus, backfill, rock,
    # both above and below insulation depth, plus the rock region below the pipe.
    gridx = [-1.0, -0.05, 0.0, 0.04, 0.08, 0.5]
    gridy = [-0.5, 0.0, 0.07]
    gridz = [5.0, 25.0, 60.0, 99.0, 101.0, 150.0]

    for i in eachindex(gridx), j in eachindex(gridy), k in eachindex(gridz)
        x, y, z = gridx[i], gridy[j], gridz[k]
        @test GeothermalWells.lookup_thermal_conductivity(acc, i, j, k, gridx, gridy, gridz) ==
              GeothermalWells.get_thermal_conductivity(x, y, z, boreholes, materials)
        @test GeothermalWells.lookup_volumetric_heat_capacity(acc, i, j, k, gridx, gridy, gridz) ==
              GeothermalWells.get_volumetric_heat_capacity(x, y, z, boreholes, materials)
    end
end

@testitem "create_cache - :auto picks accessor based on N_bh" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    backend = CPU()

    materials = HomogenousMaterialProperties{Float64}(
        2.88, 2.17e6, 0.6, 4.18e6, 44.5, 3.73e6, 0.26, 1.96e6, 1.0, 1.0
    )
    bh1 = Borehole{Float64}(0.0, 0.0, 100.0, 0.0381, 0.01, 0.0889, 0.01, 0.0989, 10.0, 50.0)
    bh2 = Borehole{Float64}(30.0, 0.0, 100.0, 0.0381, 0.01, 0.0889, 0.01, 0.0989, 10.0, 50.0)
    inlet = ConstantInlet{Float64}(20.0)

    gridx = create_adaptive_grid_1d(xmin=-60.0, xmax=60.0, dx_fine=0.01,
        growth_factor=1.5, dx_max=10.0, boreholes=(bh1, bh2),
        backend=backend, Float_used=Float64, direction=:x)
    gridy = create_adaptive_grid_1d(xmin=-60.0, xmax=60.0, dx_fine=0.01,
        growth_factor=1.5, dx_max=10.0, boreholes=(bh1, bh2),
        backend=backend, Float_used=Float64, direction=:y)
    gridz = create_uniform_gridz_with_borehole_depths(zmin=0.0, zmax=120.0, dz=10.0,
        boreholes=(bh1,), backend=backend)

    # Single borehole -> on the fly
    cache1 = create_cache(backend=backend, gridx=gridx, gridy=gridy, gridz=gridz,
        materials=materials, boreholes=(bh1,), inlet_model=inlet,
        precompute_materials=:auto)
    @test cache1.material_accessor isa GeothermalWells.OnTheFlyMaterialAccessor

    # Multiple boreholes -> precomputed
    cache2 = create_cache(backend=backend, gridx=gridx, gridy=gridy, gridz=gridz,
        materials=materials, boreholes=(bh1, bh2), inlet_model=inlet,
        precompute_materials=:auto)
    @test cache2.material_accessor isa GeothermalWells.PrecomputedMaterialAccessor
end

@testitem "create_cache - explicit precompute_materials overrides :auto" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    backend = CPU()

    materials = HomogenousMaterialProperties{Float64}(
        2.88, 2.17e6, 0.6, 4.18e6, 44.5, 3.73e6, 0.26, 1.96e6, 1.0, 1.0
    )
    bh = Borehole{Float64}(0.0, 0.0, 100.0, 0.0381, 0.01, 0.0889, 0.01, 0.0989, 10.0, 50.0)
    boreholes = (bh,)
    inlet = ConstantInlet{Float64}(20.0)

    gridx = create_adaptive_grid_1d(xmin=-60.0, xmax=60.0, dx_fine=0.01,
        growth_factor=1.5, dx_max=10.0, boreholes=boreholes,
        backend=backend, Float_used=Float64, direction=:x)
    gridy = create_adaptive_grid_1d(xmin=-60.0, xmax=60.0, dx_fine=0.01,
        growth_factor=1.5, dx_max=10.0, boreholes=boreholes,
        backend=backend, Float_used=Float64, direction=:y)
    gridz = create_uniform_gridz_with_borehole_depths(zmin=0.0, zmax=120.0, dz=10.0,
        boreholes=boreholes, backend=backend)

    base_kwargs = (; backend, gridx, gridy, gridz, materials, boreholes, inlet_model=inlet)

    # Single borehole would default to on-the-fly under :auto. Force precomputed.
    for v in (true, :precomputed)
        c = create_cache(; base_kwargs..., precompute_materials=v)
        @test c.material_accessor isa GeothermalWells.PrecomputedMaterialAccessor
    end

    for v in (false, :on_the_fly, :onthefly)
        c = create_cache(; base_kwargs..., precompute_materials=v)
        @test c.material_accessor isa GeothermalWells.OnTheFlyMaterialAccessor
    end

    @test_throws ArgumentError create_cache(; base_kwargs..., precompute_materials=:nope)
end

@testitem "create_cache - precomputed arrays equal on-the-fly evaluation" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    backend = CPU()

    materials = HomogenousMaterialProperties{Float64}(
        2.88, 2.17e6, 0.6, 4.18e6, 44.5, 3.73e6, 0.26, 1.96e6, 1.0, 1.0
    )
    bh = Borehole{Float64}(0.0, 0.0, 100.0, 0.0381, 0.01, 0.0889, 0.01, 0.0989, 10.0, 50.0)
    boreholes = (bh,)
    inlet = ConstantInlet{Float64}(20.0)

    gridx = create_adaptive_grid_1d(xmin=-30.0, xmax=30.0, dx_fine=0.01,
        growth_factor=1.5, dx_max=10.0, boreholes=boreholes,
        backend=backend, Float_used=Float64, direction=:x)
    gridy = create_adaptive_grid_1d(xmin=-30.0, xmax=30.0, dx_fine=0.01,
        growth_factor=1.5, dx_max=10.0, boreholes=boreholes,
        backend=backend, Float_used=Float64, direction=:y)
    gridz = create_uniform_gridz_with_borehole_depths(zmin=0.0, zmax=120.0, dz=10.0,
        boreholes=boreholes, backend=backend)

    cache = create_cache(backend=backend, gridx=gridx, gridy=gridy, gridz=gridz,
        materials=materials, boreholes=boreholes, inlet_model=inlet,
        precompute_materials=:precomputed)

    acc = cache.material_accessor
    @test acc isa GeothermalWells.PrecomputedMaterialAccessor
    @test size(acc.thermal_conductivity) == (length(gridz), length(gridy), length(gridx))
    @test size(acc.volumetric_heat_capacity) == (length(gridz), length(gridy), length(gridx))

    # Every entry should equal what the on-the-fly lookup computes at the same (x, y, z).
    for i in eachindex(gridx), j in eachindex(gridy), k in eachindex(gridz)
        x, y, z = gridx[i], gridy[j], gridz[k]
        @test acc.thermal_conductivity[k, j, i] ==
              GeothermalWells.get_thermal_conductivity(x, y, z, boreholes, materials)
        @test acc.volumetric_heat_capacity[k, j, i] ==
              GeothermalWells.get_volumetric_heat_capacity(x, y, z, boreholes, materials)
    end
end

@testitem "Solver: precomputed and on-the-fly produce identical results" begin
    using GeothermalWells
    using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
    using KernelAbstractions: CPU

    backend = CPU()

    materials = HomogenousMaterialProperties{Float64}(
        2.88, 2.17e6, 0.6, 4.18e6, 44.5, 3.73e6, 0.26, 1.96e6, 1.0, 1.0
    )
    bh = Borehole{Float64}(0.0, 0.0, 350.0, 0.0381, 0.01, 0.0889, 0.01, 0.0989, 10.0, 100.0)
    boreholes = (bh,)
    inlet_model = ConstantInlet{Float64}(20.0)

    gridx = create_adaptive_grid_1d(xmin=-100.0, xmax=100.0, dx_fine=0.0025,
        growth_factor=1.3, dx_max=10.0, boreholes=boreholes,
        backend=backend, Float_used=Float64, direction=:x)
    gridy = create_adaptive_grid_1d(xmin=-100.0, xmax=100.0, dx_fine=0.0025,
        growth_factor=1.3, dx_max=10.0, boreholes=boreholes,
        backend=backend, Float_used=Float64, direction=:y)
    gridz = create_uniform_gridz_with_borehole_depths(zmin=0.0, zmax=370.0, dz=10.0,
        boreholes=boreholes, backend=backend)

    T0 = initial_condition_thermal_gradient(backend, Float64, gridx, gridy, gridz;
        T_surface=2.29, gradient=0.35)

    Δt = 80.0
    tspan = (0.0, 400.0)
    saveat = [0.0, 400.0]

    function run_with(mode)
        cache = create_cache(backend=backend, gridx=gridx, gridy=gridy, gridz=gridz,
            materials=materials, boreholes=boreholes, inlet_model=inlet_model,
            precompute_materials=mode)
        prob = ODEProblem(rhs_diffusion_z!, copy(T0), tspan, cache)
        cb, sv = get_simulation_callback(saveat=saveat, print_every_n=10^9)
        solve(prob, ROCK2(max_stages=100, eigen_est=eigen_estimator),
            save_everystep=false, callback=cb, adaptive=false, dt=Δt, maxiters=Int(1e10))
        return sv, cache.T_outlet[1]
    end

    sv_pre, T_outlet_pre = run_with(:precomputed)
    sv_otf, T_outlet_otf = run_with(:on_the_fly)

    @test length(sv_pre.saveval) == 2
    @test length(sv_otf.saveval) == 2

    # The two modes only differ in whether material values are read from a buffer
    # or recomputed; they should agree to within tight floating-point tolerance.
    @test sv_pre.saveval[1] ≈ sv_otf.saveval[1] rtol=1e-10
    @test sv_pre.saveval[end] ≈ sv_otf.saveval[end] rtol=1e-8
    @test T_outlet_pre ≈ T_outlet_otf rtol=1e-8
end
