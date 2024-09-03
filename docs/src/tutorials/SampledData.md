# Clocks and Sampled-Data Systems

A sampled-data system contains both continuous-time and discrete-time components, such as a continuous-time plant model and a discrete-time control system. ModelingToolkit supports the modeling and simulation of sampled-data systems by means of *clocks*.

!!! danger "Experimental"
    
    The sampled-data interface in ModelingToolkit is currently experimental and at any time subject to breaking changes **not** respecting semantic versioning.

!!! note "Negative shifts"
    
    The initial release of the sampled-data interface **only supports negative shifts**.


## Clocks, operators and difference equations
A clock can be seen as an *event source*, i.e., when the clock ticks, an event is generated. In response to the event the discrete-time logic is executed, for example, a control signal is computed. For basic modeling of sampled-data systems, the user does not have to interact with clocks explicitly, instead, the modeling is performed using the operators

- [`ModelingToolkit.Sample`](@ref)
- [`ModelingToolkit.Hold`](@ref)
- [`ModelingToolkit.ShiftIndex`](@ref)

  or their corresponding components

- [`ModelingToolkitSampledData.Sampler`](@ref)
- [`ModelingToolkitSampledData.ZeroOrderHold`](@ref)

When a continuous-time variable `x` is sampled using `xd = Sample(x, dt)`, the result is a discrete-time variable `xd` that is defined and updated whenever the clock ticks. `xd` is *only defined when the clock ticks*, which it does with an interval of `dt`. If `dt` is unspecified, the tick rate of the clock associated with `xd` is inferred from the context in which `xd` appears. Any variable taking part in the same equation as `xd` is inferred to belong to the same *discrete partition* as `xd`, i.e., belonging to the same clock. A system may contain multiple different discrete-time partitions, each with a unique clock. This allows for modeling of multi-rate systems and discrete-time processes located on different computers etc.

To make a discrete-time variable available to the continuous partition, the `Hold` operator is used. `xc = Hold(xd)` creates a continuous-time variable `xc` that is updated whenever the clock associated with `xd` ticks, and holds its value constant between ticks.

The operators `Sample` and `Hold` are thus providing the interface between continuous and discrete partitions.

The `ShiftIndex` operator is used to refer to past and future values of discrete-time variables. The example below illustrates its use, implementing the discrete-time system

```math
x(k+1) = 0.5x(k) + u(k)
y(k) = x(k)
```

```@example clocks
using ModelingToolkit
using ModelingToolkit: t_nounits as t
@variables x(t) y(t) u(t)
dt = 0.1                # Sample interval
clock = Clock(dt)    # A periodic clock with tick rate dt
k = ShiftIndex(clock)

eqs = [
    x(k) ~ 0.5x(k - 1) + u(k - 1),
    y ~ x
]
```

A few things to note in this basic example:

  - The equation `x(k+1) = 0.5x(k) + u(k)` has been rewritten in terms of negative shifts since positive shifts are not yet supported.
  - `x` and `u` are automatically inferred to be discrete-time variables, since they appear in an equation with a discrete-time `ShiftIndex` `k`.
  - `y` is also automatically inferred to be a discrete-time-time variable, since it appears in an equation with another discrete-time variable `x`. `x,u,y` all belong to the same discrete-time partition, i.e., they are all updated at the same *instantaneous point in time* at which the clock ticks.
  - The equation `y ~ x` does not use any shift index, this is equivalent to `y(k) ~ x(k)`, i.e., discrete-time variables without shift index are assumed to refer to the variable at the current time step.
  - The equation `x(k) ~ 0.5x(k-1) + u(k-1)` indicates how `x` is updated, i.e., what the value of `x` will be at the *current* time step in terms of the *past* value. The output `y`, is given by the value of `x` at the *current* time step, i.e., `y(k) ~ x(k)`. If this logic was implemented in an imperative programming style, the logic would thus be

```julia
function discrete_step(x, u)
    x = 0.5x + u # x is updated to a new value, i.e., x(k) is computed
    y = x # y is assigned the current value of x, y(k) = x(k)
    return x, y # The state x now refers to x at the current time step, x(k), and y equals x, y(k) = x(k)
end
```

An alternative and *equivalent* way of writing the same system is

```@example clocks
eqs = [
    x(k + 1) ~ 0.5x(k) + u(k),
    y(k) ~ x(k)
]
```

but the use of positive time shifts is not yet supported. Instead, we *shifted all indices* by `-1` above, resulting in exactly the same difference equations. However, the next system is *not equivalent* to the previous one:

```@example clocks
eqs = [
    x(k) ~ 0.5x(k - 1) + u(k),
    y ~ x
]
```

In this last example, `u(k)` refers to the input at the new time point `k`., this system is equivalent to

```
eqs = [
    x(k+1) ~ 0.5x(k) + u(k+1),
    y(k) ~ x(k)
]
```

## Higher-order shifts

The expression `x(k-1)` refers to the value of `x` at the *previous* clock tick. Similarly, `x(k-2)` refers to the value of `x` at the clock tick before that. In general, `x(k-n)` refers to the value of `x` at the `n`th clock tick before the current one. As an example, the Z-domain transfer function

```math
H(z) = \dfrac{b_2 z^2 + b_1 z + b_0}{a_2 z^2 + a_1 z + a_0}
```

may thus be modeled as

```julia
@variables t y(t) [description = "Output"] u(t) [description = "Input"]
k = ShiftIndex(Clock(dt))
eqs = [
    a2 * y(k) + a1 * y(k - 1) + a0 * y(k - 2) ~ b2 * u(k) + b1 * u(k - 1) + b0 * u(k - 2)
]
```


## Initial conditions

The initial condition of discrete-time variables is defined using the `ShiftIndex` operator, for example

```julia
ODEProblem(model, [x(k-1) => 1.0], (0.0, 10.0))
```
Note how the initial condition for discrete-time variables refer to the past, and not to the value at ``t=0`` (at `tspan[1]` to be precise). The reason is that in order to perform the discrete-time state update ``x(k) = f(x(k-1), t_0)``, at the initial time ``t_0``, we need ``x(k-1)`` and not ``x(k)``.

If higher-order shifts are present, the corresponding initial conditions must be specified, e.g., the presence of the equation

```julia
x(k) = x(k - 1) + x(k - 2)
```

requires specification of the initial condition for both `x(k-1)` and `x(k-2)`.

## Multiple clocks

Multi-rate systems are easy to model using multiple different clocks. The following set of equations is valid, and defines *two different discrete-time partitions*, each with its own clock:

```julia
yd1 ~ Sample(dt1)(y)
ud1 ~ kp * (Sample(dt1)(r) - yd1)
yd2 ~ Sample(dt2)(y)
ud2 ~ kp * (Sample(dt2)(r) - yd2)
```

`yd1` and `ud1` belong to the same clock which ticks with an interval of `dt1`, while `yd2` and `ud2` belong to a different clock which ticks with an interval of `dt2`. The two clocks are *not synchronized*, i.e., they are not *guaranteed* to tick at the same point in time, even if one tick interval is a rational multiple of the other. Mechanisms for synchronization of clocks are not yet implemented.

## Accessing discrete-time variables in the solution
Discrete-time variables can be accessed with the syntax
```julia
sol[var]
```

## Implementing generic discrete-time components
Discrete-time components can be implemented without specification of the clock or sample interval. To do this, the operator `SampleTime()` is used, it returns the sample-time interval of the associated clock. See the library of components for usage example, or this simplified example component:
```julia
@mtkmodel DiscreteDerivative begin
    @extend u, y = siso = SISO()
    @parameters begin
        k = 1, [description = "Gain"]
    end
    @structural_parameters begin
        Ts = SampleTime()
        z = ShiftIndex()
    end
    @equations begin
        y(z) ~ k * (u(z) - u(z - 1)) / Ts # Ts will get the value of the associated clock, determined by clock inference.
    end
end
```
In this component, we use two structural parameters, one for the sample time and one for the shift index. By letting the sample time be a structural parameter with a default given by `SampleTime()`, the default behavior is to inherit the sample time based on the connection context (clock partition) of the component. The sample time can still be manually set in case there is no other point at which to infer the sample time, or if the user for some other reason want to override the sample time. The shift index is also a structural parameter with a default. This default leaves the clock unspecified, indicating that the clock should be inherited based on the connection context (clock partition). If the user provides a shift index with a clock specified, other components may inherit their clock from this component. The operator `SampleTime()` will in all cases return the sample time of the associated clock, no matter if this clock is explicitly specified by the user or inherited from the connection context.

In order to make components maximally generic, it is often advisable to avoid including operators like `Sample` and `Hold` at the inputs and outputs of a component, and instead let the user manually insert these operators where required. Larger components that model complete sampled-data systems may of course contain such operators internally. 

## Linearization of sampled-data systems
Linearization of discrete-time and sampled-data systems using the tools described at [Linearization of ModelingToolkit models](https://docs.sciml.ai/ModelingToolkit/dev/basics/Linearization/) is __not__ supported at the moment. 😦 

Linearization through [frequency-response analysis (FRA) is provided in JuliaSimControl](https://help.juliahub.com/juliasimcontrol/dev/linear_analysis/#Frequency-response-analysis). FRA amounts to simulating a system with wide-spectrum inputs and computing the linear small-signal transfer function using techniques from the field of system-identification.


## A complete example

Below, we model a simple continuous first-order system called `plant` that is controlled using a discrete-time controller `controller`. The reference signal is filtered using a discrete-time filter `filt` before being fed to the controller.

```@example clocks
using ModelingToolkit, Plots, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t
using ModelingToolkit: D_nounits as D
dt = 0.5 # Sample interval
@variables r(t)
clock = Clock(dt)
k = ShiftIndex(clock)

function plant(; name)
    @variables x(t)=1 u(t)=0 y(t)=0
    D = Differential(t)
    eqs = [D(x) ~ -x + u
           y ~ x]
    ODESystem(eqs, t; name = name)
end

function filt(; name) # Reference filter
    @variables x(t)=0 u(t)=0 y(t)=0
    a = 1 / exp(dt)
    eqs = [x(k) ~ a * x(k - 1) + (1 - a) * u(k)
           y ~ x]
    ODESystem(eqs, t, name = name)
end

function controller(kp; name)
    @variables y(t)=0 r(t)=0 ud(t)=0 yd(t)=0
    @parameters kp = kp
    eqs = [yd ~ Sample(y)
           ud ~ kp * (r - yd)]
    ODESystem(eqs, t; name = name)
end

@named f = filt()
@named c = controller(1)
@named p = plant()

connections = [r ~ (t >= 5)        # reference signal
               f.u ~ r             # reference to filter input
               f.y ~ c.r           # filtered reference to controller reference
               Hold(c.ud) ~ p.u    # controller output to plant input (zero-order-hold)
               p.y ~ c.y]          # plant output to controller feedback

@named cl = ODESystem(connections, t, systems = [f, c, p])
cl = complete(cl)
```

We can now simulate the system. JuliaSimCompiler is required to simulate hybrid continuous/discrete systems, we thus convert the system to an `JuliaSimCompiler.IRSystem` before calling `structural_simplify`
```@example clocks
using JuliaSimCompiler, Plots
ssys = structural_simplify(IRSystem(cl))

prob = ODEProblem(ssys, [
    cl.f.x(k-1) => 0,
], (0,15))

sol = solve(prob, Tsit5())
plot(sol, idxs=[cl.p.x, cl.c.ud], layout=(2,1))
```

Note how the initial condition provided above refers to the value of `f.x` at a past time point. If we had not provided this initial condition, we would have gotten an error like this
```@example clocks
try
    ODEProblem(ssys, [], (0,15))
catch err
    return err
end
```
The error message indicates that the -1 shift for `f.x` is not provided. The number `0.5` appearing in the error message is the period of the clock associated the variable for which an initial condition is missing `f.x`.


## Further reading
- [Analysis of linear sampled-data systems using ControlSystems.jl](https://juliacontrol.github.io/ControlSystems.jl/stable/examples/zoh/)