# Discrete block library

## Index
```@index
```

## Basic blocks
- [`Delay`](@ref)
- [`Difference`](@ref)
- [`DiscreteDerivative`](@ref)
- [`DiscreteIntegrator`](@ref)
- [`DiscreteSlewRateLimiter`](@ref)
- [`Sampler`](@ref)
- [`ZeroOrderHold`](@ref)

## Noise and measurement corruption
- [`NormalNoise`](@ref)
- [`Quantization`](@ref)
- [`SampleWithADEffects`](@ref)
- [`UniformNoise`](@ref)

## Controllers
- [`DiscretePIDParallel`](@ref)
- [`DiscretePIDStandard`](@ref)
- [`DiscreteOnOffController`](@ref)

## Discrete-time filters
- [`ExponentialFilter`](@ref)
- [`MovingAverageFilter`](@ref)


## Docstrings

```@autodocs
Modules = [ModelingToolkitSampledData]
Pages   = ["discrete_blocks.jl"]
```