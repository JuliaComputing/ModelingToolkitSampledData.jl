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