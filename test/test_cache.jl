using TestItems

@testitem "get_simulation_callback returns valid callback and saved_values" begin
    using GeothermalWells
    using DiffEqCallbacks: SavedValues

    callback, saved_values = get_simulation_callback(
        saveat=[0.0, 100.0],
        print_every_n=1000
    )

    @test saved_values isa SavedValues{Float64, Array{Float32, 3}}
    @test isempty(saved_values.t)
end

@testitem "Print callback branches" begin
    using GeothermalWells

    # Without checkpointing: callbacks are (ADI_and_ADV, print_cb, save_cb)
    callback, _ = get_simulation_callback(saveat=[0.0], print_every_n=1)
    print_cb = callback.discrete_callbacks[2]

    # t > 1 year → prints in years
    @test_nowarn print_cb.affect!((stats=(naccept=1,), t=3600.0 * 24 * 365 * 2.0))

    # t > 1 day → prints in days
    @test_nowarn print_cb.affect!((stats=(naccept=2,), t=3600.0 * 24 * 5.0))

    # t <= 1 day → prints in hours
    @test_nowarn print_cb.affect!((stats=(naccept=3,), t=7200.0))
end

@testitem "Checkpoint path helpers" begin
    using GeothermalWells

    cp = GeothermalWells._checkpoint_path("/tmp/data", "my_sim")
    @test cp == joinpath("/tmp/data", "checkpoint_my_sim.jld2")

    snap = GeothermalWells._snapshot_path("/tmp/data", "my_sim", 1)
    @test snap == joinpath("/tmp/data", "snapshot_my_sim_0001.jld2")

    snap42 = GeothermalWells._snapshot_path("/tmp/data", "my_sim", 42)
    @test snap42 == joinpath("/tmp/data", "snapshot_my_sim_0042.jld2")
end

@testitem "_clean_and_count_snapshots - empty inputs" begin
    using GeothermalWells

    tmpdir = mktempdir()

    # Empty directory
    @test GeothermalWells._clean_and_count_snapshots(tmpdir, "test") == 0

    # Nonexistent directory
    @test GeothermalWells._clean_and_count_snapshots(joinpath(tmpdir, "nope"), "test") == 0
end

@testitem "_clean_and_count_snapshots - removes orphaned snapshots" begin
    using GeothermalWells
    using JLD2: @save

    tmpdir = mktempdir()

    # Existing snapshots without a checkpoint are stale and should be removed.
    for i in 1:3
        path = GeothermalWells._snapshot_path(tmpdir, "test", i)
        u_save = zeros(Float32, 2, 2, 2)
        t_save = Float64(i)
        @save path u_save t_save
    end

    # A different checkpoint_id should be ignored and left in place.
    other_path = GeothermalWells._snapshot_path(tmpdir, "other", 1)
    u_save = ones(Float32, 2, 2, 2)
    t_save = 10.0
    @save other_path u_save t_save

    @test GeothermalWells._clean_and_count_snapshots(tmpdir, "test") == 0

    for i in 1:3
        @test !isfile(GeothermalWells._snapshot_path(tmpdir, "test", i))
    end
    @test isfile(other_path)
    @test GeothermalWells._load_existing_snapshots(tmpdir, "test", Float32)[3] == 0
end

@testitem "snapshot file matching avoids checkpoint_id prefix collisions" begin
    using GeothermalWells
    using JLD2: @save

    tmpdir = mktempdir()

    foo_path = GeothermalWells._snapshot_path(tmpdir, "foo", 1)
    u_save = fill(1.0f0, 2, 2, 2)
    t_save = 1.0
    @save foo_path u_save t_save

    foo_bar_path = GeothermalWells._snapshot_path(tmpdir, "foo_bar", 1)
    u_save = fill(2.0f0, 2, 2, 2)
    t_save = 2.0
    @save foo_bar_path u_save t_save

    times, arrays, count = GeothermalWells._load_existing_snapshots(tmpdir, "foo", Float32)
    @test count == 1
    @test times == [1.0]
    @test all(arrays[1] .== 1.0f0)

    @test GeothermalWells._clean_and_count_snapshots(tmpdir, "foo") == 0
    @test !isfile(foo_path)
    @test isfile(foo_bar_path)

    times, arrays, count = GeothermalWells._load_existing_snapshots(tmpdir, "foo_bar", Float32)
    @test count == 1
    @test times == [2.0]
    @test all(arrays[1] .== 2.0f0)
end

@testitem "_clean_and_count_snapshots - removes snapshots newer than checkpoint" begin
    using GeothermalWells
    using JLD2: @save

    tmpdir = mktempdir()

    cp_path = GeothermalWells._checkpoint_path(tmpdir, "restart")
    t_checkpoint = 200.0
    @save cp_path t_checkpoint

    for (i, t) in enumerate([0.0, 100.0, 200.0, 300.0])
        path = GeothermalWells._snapshot_path(tmpdir, "restart", i)
        u_save = fill(Float32(t), 2, 2, 2)
        t_save = t
        @save path u_save t_save
    end

    @test GeothermalWells._clean_and_count_snapshots(tmpdir, "restart") == 3

    for i in 1:3
        @test isfile(GeothermalWells._snapshot_path(tmpdir, "restart", i))
    end
    @test !isfile(GeothermalWells._snapshot_path(tmpdir, "restart", 4))

    times, arrays, count = GeothermalWells._load_existing_snapshots(tmpdir, "restart", Float32)
    @test count == 3
    @test times == [0.0, 100.0, 200.0]
    @test all(arrays[3] .== 200.0f0)
end

@testitem "_load_existing_snapshots" begin
    using GeothermalWells
    using JLD2: @save

    tmpdir = mktempdir()

    # Write 3 snapshots with known data
    for i in 1:3
        path = GeothermalWells._snapshot_path(tmpdir, "load_test", i)
        u_save = fill(Float32(i), 2, 2, 2)
        t_save = Float64(i * 100.0)
        @save path u_save t_save
    end

    times, arrays, count = GeothermalWells._load_existing_snapshots(tmpdir, "load_test", Float32)

    @test count == 3
    @test times == [100.0, 200.0, 300.0]
    @test all(arrays[1] .== 1.0f0)
    @test all(arrays[3] .== 3.0f0)
    @test eltype(arrays[1]) == Float32

    # Wrong ID returns empty
    t2, a2, c2 = GeothermalWells._load_existing_snapshots(tmpdir, "wrong_id", Float32)
    @test c2 == 0
    @test isempty(t2)
end

@testitem "_load_existing_snapshots - nonexistent directory" begin
    using GeothermalWells

    times, arrays, count = GeothermalWells._load_existing_snapshots(
        "/this/does/not/exist", "test", Float32)
    @test count == 0
    @test isempty(times)
end


@testitem "prepare_restart - fresh start" begin
    using GeothermalWells
    using KernelAbstractions: CPU

    tmpdir = mktempdir()

    T0 = fill(42.0, 3, 3, 3)
    tspan = (0.0, 1000.0)
    saveat = [0.0, 500.0, 1000.0]

    T0_out, tspan_out, saveat_out = prepare_restart(
        T0, tspan, saveat;
        checkpoint_dir=tmpdir,
        checkpoint_id="fresh",
        backend=CPU()
    )

    # No checkpoint file -> identity
    @test T0_out === T0
    @test tspan_out === tspan
    @test saveat_out === saveat
end

@testitem "prepare_restart - from checkpoint" begin
    using GeothermalWells
    using KernelAbstractions: CPU
    using JLD2: @save

    tmpdir = mktempdir()

    # Write a checkpoint at t=400
    cp_path = GeothermalWells._checkpoint_path(tmpdir, "restart")
    u_cpu = fill(99.0, 3, 3, 3)
    t_checkpoint = 400.0
    @save cp_path u_cpu t_checkpoint

    T0_fresh = fill(0.0, 3, 3, 3)
    tspan = (0.0, 1000.0)
    saveat = [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0]

    T0_out, tspan_out, saveat_out = prepare_restart(
        T0_fresh, tspan, saveat;
        checkpoint_dir=tmpdir,
        checkpoint_id="restart",
        backend=CPU()
    )

    # IC replaced with checkpoint data
    @test all(T0_out .== 99.0)

    # tspan adjusted
    @test tspan_out == (400.0, 1000.0)

    # saveat filtered to t > 400 (strict >)
    @test saveat_out == [600.0, 800.0, 1000.0]
end

@testitem "reload_snapshots!" begin
    using GeothermalWells
    using DiffEqCallbacks: SavedValues
    using JLD2: @save

    tmpdir = mktempdir()

    # Write snapshots in non-chronological file order to test sorting
    for (i, t) in enumerate([300.0, 100.0, 200.0])
        path = GeothermalWells._snapshot_path(tmpdir, "reload", i)
        u_save = fill(Float32(t), 2, 2, 2)
        t_save = t
        @save path u_save t_save
    end

    saved_values = SavedValues(Float64, Array{Float32, 3})

    reload_snapshots!(saved_values, tmpdir, "reload")

    # Should be sorted by time
    @test saved_values.t == [100.0, 200.0, 300.0]
    @test length(saved_values.saveval) == 3
    @test all(saved_values.saveval[1] .== 100.0f0)
    @test all(saved_values.saveval[3] .== 300.0f0)
end

@testitem "reload_snapshots! - empty directory" begin
    using GeothermalWells
    using DiffEqCallbacks: SavedValues

    tmpdir = mktempdir()

    saved_values = SavedValues(Float64, Array{Float32, 3})
    reload_snapshots!(saved_values, tmpdir, "empty")

    @test isempty(saved_values.t)
    @test isempty(saved_values.saveval)
end

@testitem "get_simulation_callback - without checkpointing" begin
    using GeothermalWells

    callback, saved_values = get_simulation_callback(
        saveat=[0.0, 100.0],
        print_every_n=1000
    )

    @test !isnothing(callback)
    @test isempty(saved_values.t)
end

@testitem "get_simulation_callback - with checkpointing creates directory" begin
    using GeothermalWells

    tmpdir = mktempdir()
    subdir = joinpath(tmpdir, "nested", "dir")

    callback, saved_values = get_simulation_callback(
        saveat=[0.0, 100.0],
        print_every_n=1000,
        checkpoint_dir=subdir,
        checkpoint_id="test",
        checkpoint_every_n=100
    )

    @test isdir(subdir)
end

@testitem "get_simulation_callback - removes snapshots when checkpoint missing" begin
    using GeothermalWells
    using JLD2: @save

    tmpdir = mktempdir()

    # Pre-create snapshots from an earlier run, but no checkpoint.
    for i in 1:2
        path = GeothermalWells._snapshot_path(tmpdir, "fresh", i)
        u_save = zeros(Float32, 2, 2, 2)
        t_save = Float64(i)
        @save path u_save t_save
    end

    callback, saved_values = get_simulation_callback(
        saveat=[0.0],
        checkpoint_dir=tmpdir,
        checkpoint_id="fresh",
        checkpoint_every_n=100
    )

    @test GeothermalWells._load_existing_snapshots(tmpdir, "fresh", Float32)[3] == 0
    for i in 1:2
        @test !isfile(GeothermalWells._snapshot_path(tmpdir, "fresh", i))
    end
end

@testitem "get_simulation_callback - cleans and continues snapshot numbering" begin
    using GeothermalWells
    using JLD2: @save

    tmpdir = mktempdir()

    cp_path = GeothermalWells._checkpoint_path(tmpdir, "cont")
    t_checkpoint = 3.0
    @save cp_path t_checkpoint

    # Three snapshots are covered by the checkpoint. The fourth is from after
    # the checkpoint and should be removed before continuing.
    for i in 1:4
        path = GeothermalWells._snapshot_path(tmpdir, "cont", i)
        u_save = zeros(Float32, 2, 2, 2)
        t_save = Float64(i)
        @save path u_save t_save
    end

    callback, saved_values = get_simulation_callback(
        saveat=[4.0],
        checkpoint_dir=tmpdir,
        checkpoint_id="cont",
        checkpoint_every_n=100
    )

    @test GeothermalWells._load_existing_snapshots(tmpdir, "cont", Float32)[3] == 3
    for i in 1:3
        @test isfile(GeothermalWells._snapshot_path(tmpdir, "cont", i))
    end
    @test !isfile(GeothermalWells._snapshot_path(tmpdir, "cont", 4))
end

@testitem "Checkpoint integration - restart produces same result" begin
    using GeothermalWells
    using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
    using KernelAbstractions: CPU

    backend = CPU()
    Float_used = Float64

    # --- Shared setup ---
    materials = HomogenousMaterialProperties{Float_used}(
        2.88, 2.17e6, 0.6, 4.18e6, 44.5, 3.73e6, 0.26, 1.96e6, 1.0, 1.0
    )

    borehole = Borehole{Float_used}(
        0.0, 0.0, 350.0, 0.0381, 0.01, 0.0889, 0.01, 0.0989, 10.0, 100.0
    )
    boreholes = (borehole,)

    gridx = create_adaptive_grid_1d(xmin=-100.0, xmax=100.0, dx_fine=0.0025,
        growth_factor=1.3, dx_max=10.0, boreholes=boreholes,
        backend=backend, Float_used=Float_used, direction=:x)
    gridy = create_adaptive_grid_1d(xmin=-100.0, xmax=100.0, dx_fine=0.0025,
        growth_factor=1.3, dx_max=10.0, boreholes=boreholes,
        backend=backend, Float_used=Float_used, direction=:y)
    gridz = create_uniform_gridz_with_borehole_depths(zmin=0.0, zmax=370.0, dz=10.0,
        boreholes=boreholes, backend=backend)

    T0_fresh = initial_condition_thermal_gradient(backend, Float_used, gridx, gridy, gridz;
        T_surface=2.29, gradient=0.35)

    inlet_model = ConstantInlet{Float_used}(20.0)

    cache = create_cache(backend=backend, gridx=gridx, gridy=gridy, gridz=gridz,
        materials=materials, boreholes=boreholes, inlet_model=inlet_model)

    Δt = 80.0
    tspan_full = (0.0, 800.0)
    saveat_full = [0.0, 400.0, 800.0]

    # =========================================================================
    # Run 1: Full uninterrupted simulation (reference)
    # =========================================================================
    dir_ref = mktempdir()

    T0_ref, tspan_ref, saveat_ref = prepare_restart(T0_fresh, tspan_full, saveat_full;
        checkpoint_dir=dir_ref, checkpoint_id="ref", backend=backend)

    prob_ref = ODEProblem(rhs_diffusion_z!, T0_ref, tspan_ref, cache)
    cb_ref, sv_ref = get_simulation_callback(saveat=saveat_ref,
        checkpoint_dir=dir_ref, checkpoint_id="ref", checkpoint_every_n=5)

    solve(prob_ref, ROCK2(max_stages=100, eigen_est=eigen_estimator),
        save_everystep=false, callback=cb_ref, adaptive=false, dt=Δt, maxiters=Int(1e10))

    reload_snapshots!(sv_ref, dir_ref, "ref")

    @test sv_ref.t == [0.0, 400.0, 800.0]
    @test length(sv_ref.saveval) == 3

    # =========================================================================
    # Run 2: Interrupted simulation (only first half)
    # =========================================================================
    dir_restart = mktempdir()

    tspan_half = (0.0, 480.0)  # stop after checkpoint at step 5 (t=400) but before t=800

    T0_r1, tspan_r1, saveat_r1 = prepare_restart(T0_fresh, tspan_half, saveat_full;
        checkpoint_dir=dir_restart, checkpoint_id="int", backend=backend)

    prob_r1 = ODEProblem(rhs_diffusion_z!, T0_r1, tspan_r1, cache)
    cb_r1, sv_r1 = get_simulation_callback(saveat=saveat_r1,
        checkpoint_dir=dir_restart, checkpoint_id="int", checkpoint_every_n=5)

    solve(prob_r1, ROCK2(max_stages=100, eigen_est=eigen_estimator),
        save_everystep=false, callback=cb_r1, adaptive=false, dt=Δt, maxiters=Int(1e10))

    # Should have snapshots for t=0 and t=400 on disk
    @test GeothermalWells._load_existing_snapshots(dir_restart, "int", Float32)[3] == 2

    # Checkpoint file should exist
    @test isfile(GeothermalWells._checkpoint_path(dir_restart, "int"))

    # =========================================================================
    # Run 3: Restart from checkpoint, complete remaining simulation
    # =========================================================================
    T0_r2, tspan_r2, saveat_r2 = prepare_restart(T0_fresh, tspan_full, saveat_full;
        checkpoint_dir=dir_restart, checkpoint_id="int", backend=backend)

    # Should have restarted (tspan adjusted)
    @test tspan_r2[1] > 0.0
    @test tspan_r2[2] == 800.0

    prob_r2 = ODEProblem(rhs_diffusion_z!, T0_r2, tspan_r2, cache)
    cb_r2, sv_r2 = get_simulation_callback(saveat=saveat_r2,
        checkpoint_dir=dir_restart, checkpoint_id="int", checkpoint_every_n=5)

    solve(prob_r2, ROCK2(max_stages=100, eigen_est=eigen_estimator),
        save_everystep=false, callback=cb_r2, adaptive=false, dt=Δt, maxiters=Int(1e10))

    reload_snapshots!(sv_r2, dir_restart, "int")

    # =========================================================================
    # Compare: restarted run should match reference
    # =========================================================================
    @test length(sv_r2.t) == 3
    @test sv_r2.t == sv_ref.t

    # The t=0 and t=400 snapshots should be identical (saved before the interruption)
    @test sv_r2.saveval[1] ≈ sv_ref.saveval[1]
    @test sv_r2.saveval[2] ≈ sv_ref.saveval[2]

    # The t=800 snapshot should be close but not necessarily bitwise identical,
    # since the restart checkpoint was at a different time than t=400
    # (checkpoint_every_n=5 -> every 400s, so checkpoint could be at t=400)
    @test sv_r2.saveval[3] ≈ sv_ref.saveval[3] atol=1e-6
end
