"""
Abstract base type for material property models.
"""
abstract type AbstractMaterialProperties{RealT<:Real} end

"""
Abstract base type for borehole geometries.
"""
abstract type AbstractBorehole{RealT<:Real} end

"""
    Borehole{T}(xc, yc, h, r_inner, t_inner, r_outer, t_outer, r_backfill, ṁ, insulation_depth)

U-shaped coaxial borehole heat exchanger geometry.

Fluid flows down the inner pipe, turns around at depth `h`, and returns up the outer annulus.
Precomputes squared radii and fluid velocities for performance.

# Fields
- `xc, yc`: horizontal center coordinates [m]
- `h`: borehole depth [m]
- `r_inner`: inner pipe radius [m]
- `t_inner`: inner pipe wall thickness [m]
- `r_outer`: outer pipe inner radius [m]
- `t_outer`: outer pipe wall thickness [m]
- `r_backfill`: cement backfill outer radius [m]
- `ṁ`: mass flow rate [kg/s]
- `insulation_depth`: depth to which inner pipe is insulated [m]
"""
struct Borehole{RealT<:Real} <: AbstractBorehole{RealT}
    xc::RealT
    yc::RealT
    h::RealT
    r_inner::RealT
    t_inner::RealT
    r_and_t_inner_sq::RealT # precomputed square
    r_outer::RealT
    r_outer_sq::RealT
    t_outer::RealT
    r_and_t_outer_sq::RealT # precomputed square
    r_backfill::RealT       # outer radius of cement backfill
    r_backfill_sq::RealT    # precomputed square
    ṁ::RealT               # mass flow rate
    v_inner::RealT
    v_outer::RealT
    insulation_depth::RealT

    function Borehole{T}(xc, yc, h, r_inner, t_inner, r_outer, t_outer, r_backfill, ṁ, insulation_depth) where {T<:Real}
        ρ_water = 998.2 # kg/m3
        r_and_t_inner_sq = (r_inner + t_inner)^2
        r_outer_sq = r_outer^2
        r_and_t_outer_sq = (r_outer + t_outer)^2
        r_backfill_sq = r_backfill^2
        A_inner = π * r_inner^2
        A_outer = π * r_outer^2 - π * r_and_t_inner_sq
        v_inner = ṁ / (A_inner * ρ_water)
        v_outer = ṁ / (A_outer * ρ_water)

        new{T}(T(xc), T(yc), T(h), T(r_inner), T(t_inner), T(r_and_t_inner_sq), T(r_outer), T(r_outer_sq), T(t_outer), T(r_and_t_outer_sq), T(r_backfill), T(r_backfill_sq), T(ṁ), T(v_inner), T(v_outer), T(insulation_depth))
    end
end


"""
    StratifiedMaterialProperties{N,T}

Thermal properties for stratified rock layers and borehole materials.

Supports `N` horizontal rock layers with varying properties by depth. Includes material
properties for water, steel, insulation, and cement backfill.

# Fields
- `k_rock_layers`: thermal conductivities for each rock layer [W/(m·K)]
- `rho_c_rock_layers`: volumetric heat capacities for each layer [J/(m³·K)]
- `layer_depths`: depth boundaries for each layer [m]
- `k_water`, `rho_c_water`: water thermal properties
- `k_steel`, `rho_c_steel`: steel pipe thermal properties
- `k_insulating`, `rho_c_insulating`: insulation thermal properties
- `k_backfill`, `rho_c_backfill`: cement backfill thermal properties
"""
struct StratifiedMaterialProperties{N,RealT<:Real} <: AbstractMaterialProperties{RealT}
    k_rock_layers::NTuple{N,RealT}
    rho_c_rock_layers::NTuple{N,RealT}
    layer_depths::NTuple{N,RealT}

    k_water::RealT
    rho_c_water::RealT
    k_steel::RealT
    rho_c_steel::RealT
    k_insulating::RealT
    rho_c_insulating::RealT
    k_backfill::RealT
    rho_c_backfill::RealT
end


"""
    get_rock_thermal_conductivity_by_depth(z, materials::StratifiedMaterialProperties) -> k

Return rock thermal conductivity [W/(m·K)] at depth `z` based on stratified layers.
"""
@inline function get_rock_thermal_conductivity_by_depth(z, materials::StratifiedMaterialProperties{N}) where {N}
    for i in 1:N
        if z <= materials.layer_depths[i]
            return materials.k_rock_layers[i]
        end
    end
    # If z is deeper than all specified layers, use the deepest layer
    return materials.k_rock_layers[N]
end

"""
    get_rock_volumetric_heat_capacity_by_depth(z, materials::StratifiedMaterialProperties) -> ρc

Return rock volumetric heat capacity [J/(m³·K)] at depth `z` based on stratified layers.
"""
@inline function get_rock_volumetric_heat_capacity_by_depth(z, materials::StratifiedMaterialProperties{N}) where {N}
    for i in 1:N
        if z <= materials.layer_depths[i]
            return materials.rho_c_rock_layers[i]
        end
    end
    # If z is deeper than all specified layers, use the deepest layer
    return materials.rho_c_rock_layers[N]
end

"""
    eigen_estimator_get_dmax(materials) -> d_max

Compute maximum thermal diffusivity [m²/s] across all materials.

Used for eigenvalue estimation in time stepping stability analysis.
"""
@inline function eigen_estimator_get_dmax(materials::StratifiedMaterialProperties{N}) where {N}
    # Find max diffusivity among all rock layers
    d_rock_max = zero(eltype(materials.k_rock_layers))
    for i in 1:N
        d_rock_i = materials.k_rock_layers[i] / materials.rho_c_rock_layers[i]
        d_rock_max = max(d_rock_max, d_rock_i)
    end

    return max(
        d_rock_max,
        materials.k_water / materials.rho_c_water,
        materials.k_steel / materials.rho_c_steel,
        materials.k_insulating / materials.rho_c_insulating,
        materials.k_backfill / materials.rho_c_backfill
    )
end


"""
    HomogenousMaterialProperties{T}

Thermal properties with uniform rock properties (no stratification).

Simpler alternative to [`StratifiedMaterialProperties`](@ref) when rock is homogeneous.
"""
struct HomogenousMaterialProperties{RealT<:Real} <: AbstractMaterialProperties{RealT}
    k_rock::RealT
    rho_c_rock::RealT
    k_water::RealT
    rho_c_water::RealT
    k_steel::RealT
    rho_c_steel::RealT
    k_insulating::RealT
    rho_c_insulating::RealT
    k_backfill::RealT
    rho_c_backfill::RealT
end

# Add dispatch for helper functions
@inline get_rock_thermal_conductivity_by_depth(z, materials::HomogenousMaterialProperties) = materials.k_rock

@inline get_rock_volumetric_heat_capacity_by_depth(z, materials::HomogenousMaterialProperties) = materials.rho_c_rock

@inline function eigen_estimator_get_dmax(materials::HomogenousMaterialProperties)
    return max(
        materials.k_rock / materials.rho_c_rock,
        materials.k_water / materials.rho_c_water,
        materials.k_steel / materials.rho_c_steel,
        materials.k_insulating / materials.rho_c_insulating,
        materials.k_backfill / materials.rho_c_backfill
    )
end





"""
    get_thermal_conductivity(x, y, z, boreholes, materials) -> k

Return thermal conductivity [W/(m·K)] at point `(x, y, z)`.

Determines the appropriate material (rock, water, steel, insulation, or backfill) based on
position relative to borehole geometry.
"""
@inline function get_thermal_conductivity(x, y, z, boreholes, materials)

    # Default: rock
    k = get_rock_thermal_conductivity_by_depth(z, materials)

    for bh in boreholes
        r_sq = (x - bh.xc)^2 + (y - bh.yc)^2
        # skip to next borehole if r_sq > (bh.r_outer + bh.t_outer)
        r_sq > bh.r_backfill_sq && continue

        if z <= bh.h
            if r_sq < bh.r_inner^2 # inside inner pipe
                k = materials.k_water
            elseif r_sq < bh.r_and_t_inner_sq # inner pipe wall
                k = materials.k_insulating
            elseif r_sq < bh.r_outer_sq # outer pipe
                k = materials.k_water
            elseif r_sq < bh.r_and_t_outer_sq # outer pipe wall
                if z <= bh.insulation_depth
                    k = materials.k_insulating
                else
                    k = materials.k_steel
                end
            else # r_sq < bh.r_backfill_sq # cement backfill
                k = materials.k_backfill
            end
            # small piece of insulation below pipe at "turnaround"
        elseif z <= bh.h + 3
            k = materials.k_steel
        end
    end

    return k

end


"""
    get_volumetric_heat_capacity(x, y, z, boreholes, materials) -> ρc

Return volumetric heat capacity [J/(m³·K)] at point `(x, y, z)`.

Determines the appropriate material based on position relative to borehole geometry.
"""
@inline function get_volumetric_heat_capacity(x, y, z, boreholes, materials)

    # Default: rock
    rho_c = get_rock_volumetric_heat_capacity_by_depth(z, materials)

    for bh in boreholes
        r_sq = (x - bh.xc)^2 + (y - bh.yc)^2
        # skip to next borehole if r_sq > (bh.r_outer + bh.t_outer)^2
        r_sq > bh.r_and_t_outer_sq && continue

        if z <= bh.h
            if r_sq < bh.r_inner^2 # inside inner pipe
                rho_c = materials.rho_c_water
            elseif r_sq < bh.r_and_t_inner_sq # inner pipe wall
                rho_c = materials.rho_c_insulating
            elseif r_sq < bh.r_outer_sq # outer pipe
                rho_c = materials.rho_c_water
            else # r_sq <= (bh.r_outer + bh.t_outer)^2 # outer pipe wall
                if z <= bh.insulation_depth
                    rho_c = materials.rho_c_insulating
                else
                    rho_c = materials.rho_c_steel
                end
            end
            # small piece of insulation below pipe at "turnaround"
        elseif z <= bh.h + 3
            rho_c = materials.rho_c_steel
        end
    end

    return rho_c
end





"""
Abstract base type for inlet temperature models.
"""
abstract type AbstractInletModel end

"""
    ConstantInlet{T}(T_const)

Inlet model with constant temperature.
"""
struct ConstantInlet{T} <: AbstractInletModel
    T_const::T
end
(model::ConstantInlet)(bh_idx, T_outlet_values, t) = model.T_const

"""
    HeatExchangerInlet{T}(ΔT)

Inlet model applying temperature change `ΔT` [K] to outlet temperature.

Models a heat exchanger where inlet is `T_outlet - ΔT`. Typically computed as `ΔT = Q / (ṁ * c_water)`
where `Q` [W] is heat extraction rate, `ṁ` [kg/s] is mass flow rate, and `c_water` [J/(kg·K)] is specific heat.
"""
struct HeatExchangerInlet{T} <: AbstractInletModel
    ΔT::T  # cooling/heating applied
end
(model::HeatExchangerInlet)(bh_idx, T_outlet_values, t) = T_outlet_values[bh_idx] - model.ΔT

"""
    CustomInlet{F}(func)

Inlet model with custom function `func(bh_idx, T_outlet_values, t)`.

Allows arbitrary time-varying or borehole-specific inlet temperatures.
"""
struct CustomInlet{F} <: AbstractInletModel
    func::F
end
(model::CustomInlet)(bh_idx, T_outlet_values, t) = model.func(bh_idx, T_outlet_values, t)






"""
    initial_condition_thermal_gradient(backend, Float_used, gridx, gridy, gridz; T_surface, gradient)

Create initial temperature field with linear thermal gradient.

Returns 3D array with `T(z) = T_surface + gradient * z` where `T_surface` is surface
temperature [°C] and `gradient` is thermal gradient [K/m or °C/m].
"""
function initial_condition_thermal_gradient(backend, Float_used, gridx, gridy, gridz;
                                           T_surface, gradient)
    Nx, Ny, Nz = length(gridx), length(gridy), length(gridz)
    ϕ = zeros(backend, Float_used, Nz, Ny, Nx)
    
    kernel_thermal_gradient!(backend)(ϕ, gridz, T_surface, gradient,
                                     ndrange=(Nz, Ny, Nx))
    
    return ϕ
end

@kernel function kernel_thermal_gradient!(ϕ, @Const(gridz), T_surface, gradient)
    k, j, i = @index(Global, NTuple)
    z = gridz[k]
    ϕ[k, j, i] = T_surface + gradient * z
end

# # Should be in example / eg not define in module??!?
# # Rock temperature as a function of depth [°C]
# @inline ϕ0(z::Float32) = 0.035f0 * z + 2.29f0
# @inline ϕ0(z) = 0.035 * z + 2.29

# function initial_condition_ϕ(backend, Float_used, gridx, gridy, gridz, boreholes,
#     inlet_model::AbstractInletModel, T_outlet_initial)

#     Nx, Ny, Nz = length(gridx), length(gridy), length(gridz)

#     ϕ = zeros(backend, Float_used, Nz, Ny, Nx)

#     kernel_initial_condition_ϕ!(backend)(ϕ, gridx, gridy, gridz, boreholes,
#         inlet_model, T_outlet_initial,
#         ndrange=(Nz, Ny, Nx))

#     return ϕ
# end

# @kernel function kernel_initial_condition_ϕ!(ϕ, @Const(gridx), @Const(gridy), @Const(gridz),
#     boreholes, inlet_model, @Const(T_outlet))
#     k, j, i = @index(Global, NTuple)
#     x = gridx[i]
#     y = gridy[j]
#     z = gridz[k]

#     # Default: rock temperature
#     T = ϕ0(z)

#     for (bh_idx, bh) in enumerate(boreholes)
#         r_sq = (x - bh.xc)^2 + (y - bh.yc)^2
#         r_sq > (bh.r_outer + bh.t_outer)^2 && continue

#         if z <= bh.h + 3  # includes both main section and turnaround
#             T = inlet_model(bh_idx, T_outlet, 0) #= t=0 for initial condition =#
#             break
#         end
#     end

#     ϕ[k, j, i] = T
# end