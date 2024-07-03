module ModelingToolkitSampledData
using ModelingToolkit
using JuliaSimCompiler

export get_clock
export DiscreteIntegrator, DiscreteDerivative, Delay, Difference, ZeroOrderHold, Sampler,
       ClockChanger,
       DiscretePIDParallel, DiscretePIDStandard, DiscreteStateSpace,
       DiscreteTransferFunction, NormalNoise, UniformNoise
include("discrete_blocks.jl")

end
