using TestItems

@testitem "Single borehole integration test" begin
    using GeothermalWells

    # Run the example
    example_path = joinpath(examples_dir(), "single_borehole_basic.jl")
    include(example_path)

    # Test 1: Simulation completed with expected number of snapshots
    @test length(saved_values.t) == n_saves
    @test saved_values.t[end] == tspan[2]

    # Test 2: Outlet temperature changed from initial (heat extraction occurred)
    T_outlet_final = cache.T_outlet[1]
    T_inlet = 20.0
    @test isapprox(T_outlet_final, 37.13620046384524, rtol = 0.05)

    # Test 3: Outlet temperature is physically reasonable
    T_bottom_initial = T0[end, 1, 1]
    @test T_outlet_final < T_bottom_initial  # can't exceed rock temperature
    @test T_outlet_final > T_inlet

    # Test 4: Temperature field is physically reasonable
    T_final = saved_values.saveval[end]
    @test all(isfinite, T_final)           # no NaN or Inf
    @test minimum(T_final) > -50           # no unphysical cold spots
    @test maximum(T_final) < 200           # no unphysical hot spots

    # Test 5: Surface boundary stayed roughly the same 
    @test isapprox(T_final[:, 1, 1], T0[:, 1, 1], rtol = 1e-5)
    @test isapprox(T_final[:, 1, end], T0[:, 1, end], rtol = 1e-5)
    @test isapprox(T_final[:, end, 1], T0[:, end, 1], rtol = 1e-5)
    @test isapprox(T_final[:, end, end], T0[:, end, end], rtol = 1e-5)
end