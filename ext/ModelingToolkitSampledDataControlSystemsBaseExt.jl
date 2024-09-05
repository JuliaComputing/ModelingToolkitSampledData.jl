module ModelingToolkitSampledDataControlSystemsBaseExt
using ModelingToolkitSampledData
using ControlSystemsBase: TransferFunction, issiso, numvec, denvec, Discrete


"""
    ModelingToolkitSampledData.DiscreteTransferFunction(G::TransferFunction{<:Discrete}; kwargs...)

Create a DiscreteTransferFunction from a `ControlSystems.TransferFunction`.

The sample time of `G` is used to create a periodic `Clock` object with which the transfer funciton is associated. If this is not desired, pass a custom `ShiftIndex` via the `z` keyword argument. To let the transfer function inherit the sample time of of the surrounding context (clock inference), pass and empty `z = ShiftIndex()`.
"""
function ModelingToolkitSampledData.DiscreteTransferFunction(G::TransferFunction{<:Discrete}; z = nothing, kwargs...)
    issiso(G) || throw(ArgumentError("Only SISO systems are supported"))
    b,a = numvec(G)[], denvec(G)[]
    if z === nothing
        z = ShiftIndex(Clock(G.Ts))
    end
    return DiscreteTransferFunction(b, a; z, kwargs...)
end

end