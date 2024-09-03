module ModelingToolkitSampledData
using ModelingToolkit
using JuliaSimCompiler
using StableRNGs

export get_clock
export DiscreteIntegrator, DiscreteDerivative, Delay, Difference, ZeroOrderHold, Sampler,
       ClockChanger,
       DiscretePIDParallel, DiscretePIDStandard, DiscreteStateSpace,
       DiscreteTransferFunction, NormalNoise, UniformNoise, Quantization
include("discrete_blocks.jl")

end
