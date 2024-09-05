module ModelingToolkitSampledData
using ModelingToolkit
using JuliaSimCompiler
using StableRNGs

export get_clock
export DiscreteIntegrator, DiscreteDerivative, Delay, Difference, ZeroOrderHold, Sampler,
       ClockChanger,
       DiscretePIDParallel, DiscretePIDStandard, DiscreteStateSpace,
       DiscreteTransferFunction, NormalNoise, UniformNoise, Quantization,
       ExponentialFilter
export DiscreteOnOffController
include("discrete_blocks.jl")


runtime(ssys::JuliaSimCompiler.ScheduledSystem) = ssys.discrete_info.callback.affect!.affect!.x

end
