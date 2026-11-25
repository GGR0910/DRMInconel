using Test
using DRMInconel

@testset "DRMInconel" begin
    include("test_units.jl")
    include("test_gas_properties.jl")
    include("test_equilibrium.jl")
    include("test_inconel600.jl")
    include("test_phase2.jl")
    include("test_termica.jl")
end