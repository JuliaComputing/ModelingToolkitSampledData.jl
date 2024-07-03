using ModelingToolkitSampledData
using ModelingToolkit
using JuliaSimCompiler
using Test

@testset "ModelingToolkitSampledData.jl" begin
    @testset "discrete_blocks" begin
        @info "Testing discrete_blocks"
        include("test_discrete_blocks.jl")
    end
end
