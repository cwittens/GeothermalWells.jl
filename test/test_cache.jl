using TestItems

@testitem "get_simulation_callback errors when write_to_jld=true but no data_folder_dir" begin
    using GeothermalWells

    @test_throws ErrorException get_simulation_callback(
        saveat=[0.0],
        write_to_jld=true,
        data_folder_dir=""
    )
end

@testitem "save_and_print_callback with write_to_jld=true" begin
    using GeothermalWells
    using DiffEqCallbacks
    using DiffEqBase: DiscreteCallback

    tmpdir = mktempdir()
    jld_dir = joinpath(tmpdir, "checkpoints")

    save_cb, print_cb, saved_values = GeothermalWells.save_and_print_callback(
        [0.0, 1.0];
        write_to_jld=true,
        data_folder_dir=jld_dir,
        Float_used_to_save=Float32
    )

    # Directory should have been created
    @test isdir(jld_dir)

    # Returned types
    @test print_cb isa DiscreteCallback
    @test saved_values isa DiffEqCallbacks.SavedValues{Float64, Array{Float32, 3}}
end

@testitem "save_and_print_callback with write_to_jld=false" begin
    using GeothermalWells
    using DiffEqCallbacks
    using DiffEqBase: DiscreteCallback

    save_cb, print_cb, saved_values = GeothermalWells.save_and_print_callback(
        [0.0, 1.0];
        write_to_jld=false,
        Float_used_to_save=Float32
    )

    @test print_cb isa DiscreteCallback
    @test saved_values isa DiffEqCallbacks.SavedValues{Float64, Array{Float32, 3}}
end

@testitem "print_affect! branches" begin
    using GeothermalWells

    save_cb, print_cb, saved_values = GeothermalWells.save_and_print_callback(
        [0.0];
        print_every_n=1
    )

    # Mock integrator as a NamedTuple — print_affect! reads integrator.stats.naccept and integrator.t
    # Branch: t > 3600*24*365 → prints in years
    mock_years = (stats=(naccept=1000,), t=3600.0 * 24 * 365 * 2.0)
    output_years = @test_nowarn print_cb.affect!(mock_years)

    # Branch: t > 3600*24 → prints in days
    mock_days = (stats=(naccept=2000,), t=3600.0 * 24 * 5.0)
    output_days = @test_nowarn print_cb.affect!(mock_days)

    # Branch: else (t <= 3600*24) → prints in hours
    mock_hours = (stats=(naccept=3000,), t=7200.0)
    output_hours = @test_nowarn print_cb.affect!(mock_hours)
end
