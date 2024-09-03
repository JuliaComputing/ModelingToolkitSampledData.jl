```@meta
CurrentModule = ModelingToolkitSampledData
```

# ModelingToolkitSampledData

Documentation for [ModelingToolkitSampledData](https://github.com/JuliaComputing/ModelingToolkitSampledData.jl).

ModelingToolkitSampledData.jl contains standard-library components for discrete-time systems, which when used in a model together with continuous time differential equations forms a [_sampled-data system_](https://en.wikipedia.org/wiki/Sampled_data_system).

Common examples of sampled-data systems are
- Control systems, such as when a physical device is controlled by a computer.
- Signal-processing systems.

To learn about the fundamentals of sampled-data modeling in ModelingToolkit, check out the tutorial [Clocks and Sampled-Data Systems](@ref), to see a practical example using standard-library components, have a look at [DC Motor with PI-controller](@ref). To read the docstrings of the available components, see [Discrete block library](@ref).

## Installation
To install this library, first follow the [installation instructions for JuliaSimCompiler](https://juliacomputing.github.io/JuliaSimCompiler.jl/stable/#Installing-and-Using-JuliaSimCompiler). In particular, you need to [add the JuliaHub Registry](https://help.juliahub.com/juliasim/dev/gettingstarted/juliahubregistry/). 

After the registry is added and JuliaSimCompiler is installed, you may install this package. This package is also registered in the JuliaHubRegistry, so you may add it with
```
pkg> add ModelingToolkitSampledData
```
after you have followed the steps above.
