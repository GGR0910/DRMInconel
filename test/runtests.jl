using Test
using DRMInconel

@testset "DRMInconel" begin
    include("test_units.jl")
    include("test_gas_properties.jl")
    include("test_equilibrium.jl")
end