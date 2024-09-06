using Random
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t, D_nounits as D

z = ShiftIndex()


function get_clock(s)
    ModelingToolkit.get_gui_metadata(s).type == GlobalRef(@__MODULE__, :Sampler) ||
        error("clock(sys) is supposed to be called on a ModelingToolkitStandardLibrary.Blocks.Sampler system, got $s")
    for eq in equations(s)
        td = ModelingToolkit.get_time_domain(eq.rhs)
        if td !== nothing #&& td != ModelingToolkit.AbstractClock
            return td
        end
    end
    error("No clock found")
end

"""
    DiscreteIntegrator(;name, k = 1, x = 0.0, method = :backward)

Outputs `y = ∫k*u dt`, discretized to correspond to the one of the discrete-time transfer functions
- `method = :forward`: ``T_s / (z - 1)``
- `method = :backward`: ``T_s z / (z - 1)``
- `method = :trapezoidal`: ``(T_s / 2) (z + 1) / (z - 1)``

where ``T_s`` is the sample time of the integrator.

Initial value of integrator state ``x`` can be set with `x`

# Connectors:
- `input`
- `output`

# Parameters:
 - `k`: Gain of integrator
"""
@mtkmodel DiscreteIntegrator begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        method = :backward
    end
    @variables begin
        x(t) = 0.0, [description = "State of Integrator"]
    end
    @parameters begin
        k = 1, [description = "Gain"]
    end
    @structural_parameters begin
        Ts = SampleTime()
        z = ShiftIndex()
    end
    @equations begin
        if method === :forward
            x(z) ~ x(z - 1) + k * Ts * u(z - 1)
        elseif method === :backward
            x(z) ~ x(z - 1) + k * Ts * u(z)
        elseif method ∈ (:trapezoidal, :tustin)
            x(z) ~ x(z - 1) + k * Ts * (u(z) + u(z - 1)) / 2
        end
        y ~ x(z)
    end
end

"""
    DiscreteDerivative(; k = 1)

A discrete-time derivative filter, corresponding to the transfer function
``k (z-1) / (T_s z)``

where ``T_s`` is the sample time of the derivative filter.

# Connectors:
- `input`
- `output`

# Parameters:
 - `k`: Gain of derivative filter
"""
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
        y(z) ~ k * (u(z) - u(z - 1)) / Ts
    end
end

"""
    Delay(; n = 1)

A discrete-time delay of `n` samples, corresponding to the transfer function ``z^{-n}``.

# Connectors:
- `input`
- `output`

# Structural Parameters:
 - `n`: Number of delay samples
"""
@mtkmodel Delay begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        n = 1
        z = ShiftIndex()
    end
    @equations begin
        y(z) ~ u(z - n)
    end
end

"""
    Difference()

A discrete-time difference operator, corresponding to the transfer function ``1 - z^{-1}``.

# Connectors:
- `input`
- `output`

For a discrete-time finite-difference derivative approximation, see [`DiscreteDerivative`](@ref).
"""
@mtkmodel Difference begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        z = ShiftIndex()
    end
    @equations begin
        y(z) ~ u(z) - u(z - 1)
    end
end

"""
    NormalNoise()

A discrete-time noise source that returns a normally-distributed value at each clock tick.

# Parameters
- `mu = 0`: Mean
- `sigma = 1`: Standard deviation
- `seed = 1`: Seed for the random number generator

# Structural parameters
- `z`: The `ShiftIndex` used to indicate clock partition.

# Connectors:
- `output`
"""
@mtkmodel NormalNoise begin
    @structural_parameters begin
        z = ShiftIndex()
        additive = false
    end
    @components begin
        output = RealOutput()
        if additive
            input = RealInput()
        end
    end
    @variables begin
        y(t), [description = "Output variable"]
        if additive
            u(t), [description = "Input variable"]
        end
        n(t), [description = "Internal noise variable"]
    end
    @parameters begin
        mu = 0
        sigma = 1
        seed = 1
    end
    @equations begin
        output.u ~ y
        n(z) ~ mu + sigma*Symbolics.term(seeded_randn, seed, t; type=Real)
        # n(z) ~ mu + sigma*Symbolics.term(randn, rng; type=Real)
        if additive
            y(z) ~ u(z) + n(z) + 1e-100*y(z-1) # The 0*y(z-1) is a workaround for a bug in the compiler, to force the y variable to be a discrete-time state variable
            u ~ input.u
        else
            y(z) ~ n(z) + 1e-100*y(z-1)
        end
    end
end

"""
    seeded_randn(seed, t)

Internal function. This function seeds the seed parameter as well as the current simulation time.
"""
function seeded_randn(seed, t)
    if rand() < 0.1
        @show t
    end
    rng = StableRNGs.StableRNG(hash(t, hash(seed)))
    randn(rng)
end
"""
    seeded_rand(seed, t)

Internal function. This function seeds the seed parameter as well as the current simulation time.
"""
function seeded_rand(seed, t)
    rng = StableRNGs.StableRNG(hash(t, hash(seed)))
    rand(rng)
end

"""
    UniformNoise()

A discrete-time noise source that returns a uniformly distributed value at each clock tick.

# Parameters
- `l = 0`: Lower bound
- `u = 1`: Upper bound
- `seed = 1`: Seed for the random number generator

# Structural parameters
- `rng`: A random number generator, defaults to `Random.Xoshiro()`.
- `z`: The `ShiftIndex` used to indicate clock partition.

# Connectors:
- `output`
"""
@mtkmodel UniformNoise begin
    @structural_parameters begin
        z = ShiftIndex()
        rng = Random.Xoshiro()
        additive = false
    end
    @components begin
        output = RealOutput()
        if additive
            input = RealInput()
        end
    end
    @variables begin
        y(t), [description = "Output variable"]
        n(t), [description = "Internal noise variable"]
    end
    @parameters begin
        l = 0
        u = 1
        seed = 1
    end
    @equations begin
        output.u ~ y
        n(z) ~ l + (u-l)*Symbolics.term(seeded_rand, seed, t; type=Real)
        # y(z) ~ l + (u-l)*Symbolics.term(rand, rng; type=Real)

        if additive
            y(z) ~ input.u(z) + n(z) + 1e-100*y(z-1) # The 0*y(z-1) is a workaround for a bug in the compiler, to force the y variable to be a discrete-time state variable
        else
            y(z) ~ n(z) + 1e-100*y(z-1)
        end
    end
end

"""
    GenericNoise()

A discrete-time noise source that at each clock tick returns a random value distributed according to the provided distribution.

# Structural parameters
- `rng`: A random number generator, defaults to `Random.Xoshiro()`.
- `z`: The `ShiftIndex` used to indicate clock partition.
- `d`: The distribution to sample from.`

# Connectors:
- `output`
"""
# @mtkmodel GenericNoise begin
#     @components begin
#         output = RealOutput()
#     end
#     @variables begin
#         y(t), [description = "Output variable"]
#     end
#     @structural_parameters begin
#         z = ShiftIndex()
#         rng = Random.Xoshiro()
#         d
#     end
#     @equations begin
#         y(z) ~ Symbolics.term(rand, rng, d; type=Real)
#         output.u ~ y
#     end
# end

"""
    ZeroOrderHold()

Zero-order-Hold translates a discrete time (clocked) signal into continuous time by holding the value of the discrete signal constant until the next sample.

# Connectors:
- `input` (discrete-time signal)
- `output` (continuous-time signal)
"""
@mtkmodel ZeroOrderHold begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        z = ShiftIndex()
    end
    @equations begin
        y ~ Hold(u(z))
    end
end

"""
    Sampler()
    Sampler(; dt::Real)
    Sampler(; clock::AbstractClock)

`Sampler` transforms a continuous-time signal into discrete time by sampling the input signal every time the associated clock ticks. The clock can be specified explicitly using the `clock` keyword argument, or implicitly by providing a sample time `dt`, in which case a standard periodic `Clock` is used. 

# Connectors:
- `input` (continuous-time signal)
- `output` (discrete-time signal)
"""
@mtkmodel Sampler begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        dt = nothing
        clock = (dt === nothing ? InferredDiscrete() : Clock(dt))
    end
    @equations begin
        y ~ Sample(clock)(u)
    end
end

"""
    ClockChanger()
    ClockChanger(; to::AbstractClock, from::AbstractClock)

# Connectors:
- `input` (continuous-time signal)
- `output` (discrete-time signal)

!!! danger "Experimental"
    The ClockChanger component is experimental and has known correctness issues. Please use with caution.
"""
@mtkmodel ClockChanger begin
    begin
        isdefined(Main, :JuliaSimCompiler) || error("JuliaSimCompiler must be defined in the Main module for the ClockChanger component to work. Run `import JuliaSimCompiler`.")
        @warn "The ClockChanger component is experimental and has known correctness issues. Please use with caution."
    end
    @extend u, y = siso = SISO()
    @structural_parameters begin
        to
        from
    end
    @parameters begin
        O = 0 # This is just a dummy to workaround a bug in JSComp
    end
    @equations begin
        y(ShiftIndex(to)) ~ Main.JuliaSimCompiler.ClockChange(; to, from)(u(ShiftIndex(from))) + O
    end
end

"""
    DiscretePIDParallel(;name, kp = 1, ki = 1, kd = 1, Ni = √(max(kd * ki, 1e-6)), Nd = 10kp, u_max = Inf, u_min = -u_max, wp = 1, wd = 1, Ts = 1, with_I = true, with_D = true, Imethod = :forward, Dmethod = :backward)

Discrete-time PID controller on parallel form with anti-windup and set-point weighting.

The controller is implemented on parallel form:

Simplified:
```math
u = k_p e + \\int k_i e dt + k_d \\dfrac{de}{dt} 
```

Detailed:
```math
u = k_p(w_p r - y) + \\int \\big( k_i (r - y) + N_i e_s \\big ) dt + k_d \\dfrac{d}{dt}(w_d r - y)
```

where `e_s = u - v` is the saturated error signal, `v` is the unsaturated control signal and `u` is the saturated control signal.

The derivative is filtered to allow a maximum gain of ``N_d``.

The integrator is discretized using the method specified by `Imethod`, options include
- `Imethod = :forward` (default): Corresponding to the transfer function ``T_s / (z - 1)``
- `Imethod = :backward`: Corresponding to the transfer function ``T_s z / (z - 1)``
- `Imethod = :trapezoidal`: Corresponding to the transfer function ``(T_s / 2) (z + 1) / (z - 1)``

The derivative is discretized using the method specified by `Dmethod`, options include
- `Dmethod = :forward`: Corresponding to the transfer function ``\\dfrac{N (z-1)}{z - \\dfrac{k_d-N T_s}{k_d}}``.
- `Dmethod = :backward` (default): Corresponding to the transfer function ``\\dfrac{\\dfrac{Nk_d}{k_d + N T_s}(z-1)}{z - \\dfrac{k_d}{k_d + N T_s}}``

Anti windup is realized by tracking using the gain ``N_i`` on the error signal ``e_s`` when the output is saturated.

To use the controller in 1DOF mode, i.e., with only the control error as input, connect the error signal to the `reference` connector, connect a `Constant(; k = 0)` to the `measurement` connector and set `wp = wd = 1`.

# Connectors:
- `reference`: The reference signal to the controller (or the error signal if used in 1DOF mode)
- `measurement`: The measurement feedback
- `ctr_output`: The control signal output

# Parameters:
- `kp`: Proportional gain
- `ki`: Integral gain (only active if `with_I = true`)
- `kd`: Derivative gain (only active if `with_D = true`)
- `Ni`: Anti-windup gain (only active if `with_I = true`)
- `Nd`: Maximum derivative gain (only active if `with_D = true`). Typically set to 10-100 times the proportional gain.
- `u_max`: Maximum output above which the output is saturated
- `u_min`: Minimum output below which the output is saturated. This defaults to `-u_max` if `u_max > 0` and `-Inf` otherwise.
- `wp`: `[0, 1]` Set-point weighting in the proportional part. Set to `0` to prevent step changes in the output due to step changes in the reference.
- `wd`: `[0, 1]` Set-point weighting in the derivative part. Set to `0` to prevent very large impulsive changes in the output due to step changes in the reference.
- `with_I`: Whether or not to include the integral part
- `with_D`: Whether or not to include the derivative part
- `Imethod`: Discretization method for the integrator (see details above)
- `Dmethod`: Discretization method for the derivative (see details above)


# Extended help:
## Internal variables:
- `I`: State of integrator
- `D`: State of filtered derivative
- `r`: Reference signal internal variable
- `y`: Measurement signal internal variable
- `wde`: Setpoint-weighted error for derivative
- `v`: Un-saturated output of the controller
- `u`: Saturated output of the controller
- `eI`: Error signal input to integrator including anit-windup tracking signal
- `e`: Error signal
"""
@mtkmodel DiscretePIDParallel begin
    @structural_parameters begin
        Imethod = :forward
        Dmethod = :backward
        with_I = true
        with_D = true
        Ts = SampleTime()
        z = ShiftIndex()
    end
    @components begin
        reference = RealInput()
        measurement = RealInput()
        ctr_output = RealOutput()
    end
    @variables begin
        I(t) = 0.0, [description = "State of Integrator"]
        D(t) = 0.0, [description = "State of filtered derivative"]
        r(t), [guess = 0, description = "Reference signal internal variable"]
        y(t), [guess = 0, description = "Measurement signal internal variable"]
        wde(t), [guess = 0, description = "Setpoint-weighted error for derivative"]
        v(t), [guess = 0, description = "Un-saturated output of the controller"]
        u(t), [guess = 0, description = "Saturated output of the controller"]
        eI(t),
        [guess = 0,
            description = "Error signal input to integrator including anit-windup tracking signal"]
        e(t), [guess = 0, description = "Error signal"]
    end
    @parameters begin
        kp = 1, [description = "Proportional gain"]
        ki = 1, [description = "Integral gain"]
        kd = 1, [description = "Derivative gain"]
        Ni = √(max(kd * ki, 1e-6)), [description = "Anti-windup gain"]
        Nd = 10 * kp, [description = "Maximum derivative gain"]
        u_max = Inf, [description = "Maximum output"]
        u_min = ifelse(u_max > 0, -u_max, -Inf), [description = "Minimum output"]
        wp = 1, [description = "Set-point weighting in the proportional part."]
        wd = 1, [description = "Set-point weighting in the derivative part."]
    end
    @equations begin
        r ~ reference.u
        y ~ measurement.u
        u ~ ctr_output.u
        e ~ r - y
        v ~ kp * (wp * r - y) + I(z - 1) + D # Unsaturated control signal
        u ~ Blocks._clamp(v, u_min, u_max) # Saturated control signal
        if with_I
            eI ~ e + Ni * (u - v) # Add anti-windup tracking signal to error before integration
            if Imethod === :forward
                I(z) ~ I(z - 1) + Ts * ki * eI(z - 1)
            elseif Imethod === :backward
                I(z) ~ I(z - 1) + Ts * ki * eI(z)
            elseif Imethod ∈ (:trapezoidal, :tustin)
                I(z) ~ I(z - 1) + Ts * ki * (eI(z) + eI(z - 1)) / 2
            else
                error("Unknown integrator discretization method $Imethod, must be one of :forward, :backward, :trapezoidal")
            end
        else
            I(z) ~ 0
        end
        if with_D
            wde ~ wd * r - y
            if Dmethod === :forward
                D(z) ~ (kd - Nd * Ts) / kd * D(z - 1) + Nd * (wde(z) - wde(z - 1))
            elseif Dmethod === :backward
                D(z) ~ kd / (kd + Nd * Ts) * D(z - 1) +
                       Nd * kd / (kd + Nd * Ts) * (wde(z) - wde(z - 1))
            else
                error("Unknown derivative discretization method $Dmethod, must be one of :forward, :backward")
            end
        else
            D(z) ~ 0
        end
    end
end

"""
    DiscretePIDStandard(;name, K = 1, Ti = 1, Td = 1, Ni = √(max(kd * ki, 1e-6)), Nd = 10, u_max = Inf, u_min = -u_max, wp = 1, wd = 1, Ts = 1, with_I = true, with_D = true, Imethod = :forward, Dmethod = :backward)

Discrete-time PID controller on standard (ideal) form with anti-windup and set-point weighting.

The controller is implemented on standard form

Simplified:
```math
u = K \\left( e + \\int \\frac{1}{T_i} e  dt + T_d \\dfrac{de}{dt} \\right)
```

Detailed:
```math
u = K \\left( (w_p r - y) + \\int \\big( \\frac{1}{T_i} (r - y) + N_i e_s \\big ) dt + T_d \\dfrac{d}{dt}(w_d r - y) \\right)
```

where `e_s = u - v` is the saturated error signal, `v` is the unsaturated control signal and `u` is the saturated control signal.

The derivative is filtered to allow a maximum gain of ``N_d``.

The integrator is discretized using the method specified by `Imethod`, options include
- `Imethod = :forward` (default): Corresponding to the transfer function ``T_s / (z - 1)``
- `Imethod = :backward`: Corresponding to the transfer function ``T_s z / (z - 1)``
- `Imethod = :trapezoidal`: Corresponding to the transfer function ``(T_s / 2) (z + 1) / (z - 1)``

The derivative is discretized using the method specified by `Dmethod`, options include
- `Dmethod = :forward`: Corresponding to the transfer function ``\\dfrac{N (z-1)}{z - \\dfrac{T_d-N T_s}{T_d}}``.
- `Dmethod = :backward` (default): Corresponding to the transfer function ``\\dfrac{\\dfrac{NT_d}{T_d + N T_s}(z-1)}{z - \\dfrac{T_d}{T_d + N T_s}}``

Anti windup is realized by tracking using the gain ``N_i`` on the error signal ``e_s`` when the output is saturated.

To use the controller in 1DOF mode, i.e., with only the control error as input, connect the error signal to the `reference` connector, connect a `Constant(; k = 0)` to the `measurement` connector and set `wp = wd = 1`.

# Connectors:
- `reference`: The reference signal to the controller (or the error signal if used in 1DOF mode)
- `measurement`: The measurement feedback
- `ctr_output`: The control signal output

# Parameters:
- `K`: Proportional gain
- `Ti`: Integral time constant (only active if `with_I = true`)
- `Td`: Derivative time (only active if `with_D = true`)
- `Ni`: Anti-windup gain (only active if `with_I = true`)
- `Nd`: Maximum derivative gain (only active if `with_D = true`). Typically set to 10-100.
- `u_max`: Maximum output above which the output is saturated
- `u_min`: Minimum output below which the output is saturated. This defaults to `-u_max` if `u_max > 0` and `-Inf` otherwise.
- `wp`: `[0, 1]` Set-point weighting in the proportional part. Set to `0` to prevent step changes in the output due to step changes in the reference.
- `wd`: `[0, 1]` Set-point weighting in the derivative part. Set to `0` to prevent very large impulsive changes in the output due to step changes in the reference.
- `with_I`: Whether or not to include the integral part
- `with_D`: Whether or not to include the derivative part
- `Imethod`: Discretization method for the integrator (see details above)
- `Dmethod`: Discretization method for the derivative (see details above)


# Extended help:
## Internal variables:
- `I`: State of integrator
- `D`: State of filtered derivative
- `r`: Reference signal internal variable
- `y`: Measurement signal internal variable
- `wde`: Setpoint-weighted error for derivative
- `v`: Un-saturated output of the controller
- `u`: Saturated output of the controller
- `eI`: Error signal input to integrator including anit-windup tracking signal
- `e`: Error signal
"""
@mtkmodel DiscretePIDStandard begin
    @structural_parameters begin
        Imethod = :forward
        Dmethod = :backward
        with_I = true
        with_D = true
        Ts = SampleTime()
        z = ShiftIndex()
    end
    @components begin
        reference = RealInput()
        measurement = RealInput()
        ctr_output = RealOutput()
    end
    @variables begin
        I(t) = 0.0, [description = "State of Integrator"]
        D(t) = 0.0, [description = "State of filtered derivative"]
        r(t) = 0.0, [description = "Reference signal internal variable"]
        y(t) = 0.0, [description = "Measurement signal internal variable"]
        wde(t) = 0.0, [description = "Setpoint-weighted error for derivative"]
        v(t) = 0.0, [description = "Un-saturated output of the controller"]
        u(t) = 0.0, [description = "Saturated output of the controller"]
        eI(t) = 0.0,
        [description = "Error signal input to integrator including anit-windup tracking signal"]
        e(t) = 0.0, [description = "Error signal"]
    end
    @parameters begin
        K = 1, [description = "Proportional gain"]
        Ti = 1, [description = "Integral time"]
        Td = 1, [description = "Derivative time"]
        Ni = √(max(Td / Ti, 1e-6)), [description = "Anti-windup gain"]
        Nd = 10, [description = "Maximum derivative gain"]
        u_max = Inf, [description = "Maximum output"]
        u_min = ifelse(u_max > 0, -u_max, -Inf), [description = "Minimum output"]
        wp = 1, [description = "Set-point weighting in the proportional part."]
        wd = 1, [description = "Set-point weighting in the derivative part."]
    end
    @equations begin
        r ~ reference.u
        y ~ measurement.u
        u ~ ctr_output.u
        e ~ r - y
        v ~ K * ((wp * r - y) + I(z - 1) + D) # Unsaturated control signal
        u ~ Blocks._clamp(v, u_min, u_max) # Saturated control signal
        if with_I
            eI ~ e + Ni * (u - v) # Add anti-windup tracking signal to error before integration
            if Imethod === :forward
                I(z) ~ I(z - 1) + Ts / Ti * eI(z - 1)
            elseif Imethod === :backward
                I(z) ~ I(z - 1) + Ts / Ti * eI(z)
            elseif Imethod ∈ (:trapezoidal, :tustin)
                I(z) ~ I(z - 1) + Ts / Ti * (eI(z) + eI(z - 1)) / 2
            else
                error("Unknown integrator discretization method $Imethod, must be one of :forward, :backward, :trapezoidal")
            end
        else
            I(z) ~ 0
        end
        if with_D
            wde = wd * r - y
            if Dmethod === :forward
                D(z) ~ (Td - Nd * Ts) / Td * D(z - 1) + Nd * (wde(z) - wde(z - 1))
            elseif Dmethod === :backward
                D(z) ~ Td / (Td + Nd * Ts) * D(z - 1) +
                       Nd * Td / (Td + Nd * Ts) * (wde(z) - wde(z - 1))
            else
                error("Unknown derivative discretization method $Dmethod, must be one of :forward, :backward")
            end
        else
            D(z) ~ 0
        end
    end
end

# """
#     DiscreteStateSpace(; A, B, C, D, u0 = zeros(size(B, 2)), y0 = zeros(size(C, 1)))

# A linear, time invariant, discrete-time state-space system on the form
# ```math
# x_{k+1} = A x_k + B u_k
# y_k = C x_k + D u_k
# ```

# `y0` and `u0` can be used to set an operating point, providing them changes the dynamics from an LTI system to the affine system
# ```math
# x_{k+1} = A x_k + B (u_k - u_0)
# y_k = C x_k + D (u_k - u_0) + y_0
# ```
# """
# @mtkmodel DiscreteStateSpace begin
#     @structural_parameters begin
#         z = ShiftIndex()
#     end
#     @parameters begin
#         A, [description = "State matrix"]
#         B, [description = "Input matrix"]
#         C, [description = "Output matrix"]
#         (D = 0),      [description = "Feedthrough matrix"]
#         (u0[1:size(B, 2)] = 0), [description = "Input operating point"]
#         (y0[1:size(C, 1)] = 0), [description = "Output operating point"]
#     end
#     begin
#         nx = size(A, 1)
#         nu = size(B, 2)
#         ny = size(C, 1)
#         size(A, 2) == nx || error("`A` has to be a square matrix.")
#         size(B, 1) == nx || error("`B` has to be of dimension ($nx x $nu).")
#         size(C, 2) == nx || error("`C` has to be of dimension ($ny x $nx).")
#     end
#     @components begin
#         input = RealInput(nin = nu)
#         output = RealOutput(nout = ny)
#     end
#     @variables begin
#         (x(t)[1:nx]), [description = "State"]
#         (u(t)[1:nu]), [description = "Input"]
#         (y(t)[1:ny]), [description = "Output"]
#     end
#     @equations begin
#         x(z) ~ A * x(z-1) .+ B * (u(z-1) .- u0)
#         y(z) ~ C * x(z) .+ D * (u(z) .- u0) .+ y0
#         output.u ~ y
#         input.u ~ u
#     end
# end

# function DiscreteStateSpace(A, B, C, D = nothing; kwargs...)
#     DiscreteStateSpace(; A, B, C, D, kwargs...)
# end

"""
    DiscreteTransferFunction(; b, a)

A discrete-time transfer function on the form
```math
H(z) = \\dfrac{B(z)}{A(z)} = \\dfrac{b_{n_b} z^{n_b} + b_{n_b-1} z^{n_b-1} + \\cdots + b_1 z + b_0}{a_{n_a} z^{n_a} + a_{n_a-1} z^{n_a-1} + \\cdots + a_1 z + a_0}
```

With the coeffiencents specified in decreasing orders of ``z``, i.e., ``b = [b_{n_b}, b_{n_b-1}, \\cdots, b_1, b_0]`` and ``a = [a_{n_a}, a_{n_a-1}, \\cdots, a_1, a_0]``.

## Parameters:
- `b`: Numerator coefficients of transfer function (e.g., z - 1 is specified as `[1,-1]`)
- `a`: Denominator coefficients of transfer function (e.g., z^2 - 0.78z + 0.37 is specified as `[1, -0.78, 0.37]`)
- `verbose`: Whether or not to print a warning if the numerator degree is larger than the denominator degree.

## Connectors:
- `input`: Input signal
- `output`: Output signal

# Extended help:
This component supports SISO systems only. To simulate MIMO transfer functions, use [ControlSystemsBase.jl](https://juliacontrol.github.io/ControlSystems.jl/stable/man/creating_systems/) to convert the transfer function to a statespace system, optionally compute a minimal realization using [`minreal`](https://juliacontrol.github.io/ControlSystems.jl/stable/lib/constructors/#ControlSystemsBase.minreal), and then use a [`DiscreteStateSpace`](@ref) component instead.

See also [ControlSystemsMTK.jl](https://juliacontrol.github.io/ControlSystemsMTK.jl/stable/) for an interface between [ControlSystems.jl](https://juliacontrol.github.io/ControlSystems.jl/stable/) and ModelingToolkit.jl for advanced manipulation of transfer functions and linear statespace systems. For linearization, see [`linearize`](@ref) and [Linear Analysis](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/linear_analysis/).
"""
# @mtkmodel DiscreteTransferFunction begin
#     @parameters begin
#         b = [1], [description = "Numerator coefficients of transfer function (e.g., z - 1 is specified as [1,-1])"]
#         a = [1], [description = "Denominator coefficients of transfer function (e.g., z^2 - 0.78z + 0.37 is specified as [1, -0.78, 0.37])"]
#     end
#     @structural_parameters begin
#         verbose = true
# Ts = SampleTime()
# z = ShiftIndex()
#     end
#     begin
#         na = length(a)
#         nb = length(b)
#         verbose && nb > na && @warn "DiscreteTransferFunction: Numerator degree is larger than denominator degree, this is not a proper transfer function. Simulation of a model including this transfer funciton will require at least $(nb-na) samples additional delay in series. Silence this warning with verbose=false"
#         Ts = SampleTime()
#     end
#     @components begin
#         input = RealInput()
#         output = RealOutput()
#     end
#     @variables begin
#         u(t) = 0.0, [description = "Input flowing through connector `input`"]
#         y(t) = 0.0, [description = "Output flowing through connector `output`"]
#     end
#     @equations begin
#         sum(y(z+k-1) * a[na-k+1] for k in 1:na) ~ sum(u(z+k-1) * b[nb-k+1] for k in 1:nb)
#         input.u ~ u
#         output.u ~ y
#     end
# end

# DiscreteTransferFunction(b, a; kwargs...) = DiscreteTransferFunction(; b, a, kwargs...)

##



# https://github.com/SciML/ModelingToolkit.jl/issues/2843
# @component function DiscreteStateSpace(; A, B, C, D = nothing, x = zeros(size(A, 1)), name,
#         u0 = zeros(size(B, 2)), y0 = zeros(size(C, 1)), z = z)
#     nx, nu, ny = size(A, 1), size(B, 2), size(C, 1)
#     size(A, 2) == nx || error("`A` has to be a square matrix.")
#     size(B, 1) == nx || error("`B` has to be of dimension ($nx x $nu).")
#     size(C, 2) == nx || error("`C` has to be of dimension ($ny x $nx).")
#     if B isa AbstractVector
#         B = reshape(B, length(B), 1)
#     end
#     if isnothing(D) || iszero(D)
#         D = zeros(ny, nu)
#     else
#         size(D) == (ny, nu) || error("`D` has to be of dimension ($ny x $nu).")
#     end
#     @named input = RealInput(nin = nu)
#     @named output = RealOutput(nout = ny)
#     @variables x(t)[1:nx]=x [
#         description = "State variables"
#     ]
#     x = collect(x)
#     # pars = @parameters A=A B=B C=C D=D # This is buggy
#     eqs = [ 
#         [x[i](z) ~ sum(A[i, k] * x[k](z-1) for k in 1:nx) +
#                                  sum(B[i, j] * (input.u[j](z-1) - u0[j]) for j in 1:nu)
#          for i in 1:nx]; # cannot use D here
#         [output.u[j] ~ sum(C[j, i] * x[i] for i in 1:nx) +
#                        sum(D[j, k] * (input.u[k] - u0[k]) for k in 1:nu) + y0[j]
#          for j in 1:ny];
#     ]
#     @show eqs
#     compose(ODESystem(eqs, t, name = name), [input, output])
# end

# DiscreteStateSpace(A, B, C, D = nothing; kwargs...) = DiscreteStateSpace(; A, B, C, D, kwargs...)

# symbolic_eps(t) = eps(t)
# @register_symbolic symbolic_eps(t)

# """
#     TransferFunction(; b, a, name)

# A single input, single output, linear time-invariant system provided as a transfer-function.
# ```
# Y(s) = b(s) / a(s)  U(s)
# ```
# where `b` and `a` are vectors of coefficients of the numerator and denominator polynomials, respectively, ordered such that the coefficient of the highest power of `s` is first.

# The internal state realization is on controller canonical form, with state variable `x`, output variable `y` and input variable `u`. For numerical robustness, the realization used by the integrator is scaled by the last entry of the `a` parameter. The internally scaled state variable is available as `x_scaled`.

# To set the initial state, it's recommended to set the initial condition for `x`, and let that of `x_scaled` be computed automatically.

# # Parameters:
# - `b`: Numerator polynomial coefficients, e.g., `2s + 3` is specified as `[2, 3]`
# - `a`: Denominator polynomial coefficients, e.g., `s² + 2ωs + ω^2` is specified as `[1, 2ω, ω^2]`

# # Connectors:
#   - `input`
#   - `output`

# See also [`StateSpace`](@ref) which handles MIMO systems, as well as [ControlSystemsMTK.jl](https://juliacontrol.github.io/ControlSystemsMTK.jl/stable/) for an interface between [ControlSystems.jl](https://juliacontrol.github.io/ControlSystems.jl/stable/) and ModelingToolkit.jl for advanced manipulation of transfer functions and linear statespace systems. For linearization, see [`linearize`](@ref) and [Linear Analysis](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/linear_analysis/).
# """
# @component function DiscreteTransferFunction(; b = [1], a = [1, 1], name, z=z)
#     nb = length(b)
#     na = length(a)
#     nb <= na ||
#         error("Transfer function is not proper, the numerator must not be longer than the denominator")
#     nx = na - 1
#     nbb = max(0, na - nb)

#     @named begin
#         input = RealInput()
#         output = RealOutput()
#     end

#     @parameters begin
#         b[1:nb] = b,
#         [
#             description = "Numerator coefficients of transfer function (e.g., 2s + 3 is specified as [2,3])"
#         ]
#         a[1:na] = a,
#         [
#             description = "Denominator coefficients of transfer function (e.g., `s² + 2ωs + ω^2` is specified as [1, 2ω, ω^2])"
#         ]
#         bb[1:(nbb + nb)] = [zeros(nbb); b]
#         d = bb[1] / a[1], [description = "Direct feedthrough gain"]
#     end

#     a = collect(a)
#     @parameters a_end = ifelse(a[end] > 100 * symbolic_eps(sqrt(a' * a)), a[end], 1.0)

#     pars = [collect(b); a; collect(bb); d; a_end]
#     @variables begin
#         x(t)[1:nx] = zeros(nx),
#         [description = "State of transfer function on controller canonical form"]
#         x_scaled(t)[1:nx] = collect(x) * a_end, [description = "Scaled vector x"]
#         u(t), [description = "Input flowing through connector `input`"]
#         y(t), [description = "Output flowing through connector `output`"]
#     end

#     x = collect(x)
#     x_scaled = collect(x_scaled)
#     bb = collect(bb)

#     sts = [x; x_scaled; y; u]

#     if nx == 0
#         eqs = [y ~ d * u]
#     else
#         eqs = Equation[x_scaled[1](z) ~ (-a[2:na]'x_scaled(z-1) + a_end * u) / a[1]
#                        [x_scaled[i](z) .~ x_scaled[j](z-1) for (i,j) in zip(2:nx, 1:(nx - 1))]
#                        y ~ ((bb[2:na] - d * a[2:na])'x_scaled) / a_end + d * u
#                        x .~ x_scaled ./ a_end]
#     end
#     push!(eqs, input.u ~ u)
#     push!(eqs, output.u ~ y)
#     compose(ODESystem(eqs, t, sts, pars; name = name), input, output)
# end


"""
    DiscreteSlewRateLimiter(rate = 1.0, rate_negative = rate)

A discrete-time slew rate limiter that limits the rate of change of the input signal.

Note, the sample interval is not taken into account when computing the rate of change, the difference between two consequetive samples is saturated.

# Parameters:
- `rate`: Slew rate limit (in positive/increasing direction). Must be a positive number.
- `rate_negative`: Negative slew rate limit, defaults to `rate`. Must be a positive number.

# Variables
- `d`: Unsaturated rate of change of the input signal
- `u`: Input signal
- `y`: Output signal (saturated slew rate)

# Connectors:
- `input`
- `output`

# Example
```
cl = Clock(0.1)
z = ShiftIndex(cl)
@mtkmodel SlweRateLimiterModel begin
    @components begin
        input = Sine(amplitude=1, frequency=0.8)
        limiter = DiscreteSlewRateLimiter(; z, rate=0.4, rate_negative = 0.3)
    end
    @variables begin
        x(t) = 0 # Dummy variable to workaround JSCompiler bug
    end
    @equations begin
        connect(input.output, limiter.input)
        D(x) ~ 0.1x + Hold(limiter.y)
    end
end
@named m = SlweRateLimiterModel()
m = complete(m)
ssys = structural_simplify(IRSystem(m))
prob = ODEProblem(ssys, [m.limiter.y(z-1) => 0], (0.0, 2.0))
sol = solve(prob, Tsit5(), dtmax=0.01)
plot(sol, idxs=[m.input.output.u, m.limiter.y], title="Slew rate limited sine wave")
```
"""
@mtkmodel DiscreteSlewRateLimiter begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        z = ShiftIndex()
    end
    @parameters begin
        rate = 1.0, [description = "Slew rate limit"]
        rate_negative = rate, [description = "Negative slew rate limit"]
    end
    @variables begin
        d(t)
    end
    @equations begin
        d(z) ~ u(z) - y(z-1)
        y(z) ~ y(z-1) + clamp(d(z), -rate_negative, rate)
    end
end

"""
    Quantization

A quantization block that quantizes the input signal to a specified number of bits.

# Parameters:
- `y_max`: Upper limit of output
- `y_min`: Lower limit of output
- `bits`: Number of bits of quantization
- `quantized`: If quantization effects shall be computed. If false, the output is equal to the input, which may be useful for, e.g., linearization.
- `midrise`: (structural) If true (default), the quantizer is a midrise quantizer, otherwise it is a midtread quantizer. See [Docs: Quantization](https://juliacomputing.github.io/ModelingToolkitSampledData.jl/dev/tutorials/noise/#Quantization) for more details.

# Connectors:
- `input`
- `output`

# Variables
- `y`: Output signal, equal to `output.u`
- `u`: Input signal, equal to `input.u`

# Automatic differentiation
This block is not differentiable, the derivative is zero everywhere exect for at the level transition where it is ill defined. To use in a differentiable context, set `quantized = false` which turns this block into the identity function.
"""
@mtkmodel Quantization begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        z = ShiftIndex()
        midrise = true
    end
    @parameters begin
        y_max = 1, [description = "Upper limit of output"]
        y_min = -1, [description = "Lower limit of output"]
        bits = 8, [description = "Number of bits of quantization"]
        quantized::Bool = true, [description = "If quantization effects shall be computed."]
    end
    @equations begin
        y(z) ~ ifelse(quantized == true, quantize(u(z), bits, y_min, y_max, midrise), u(z))
    end
end

function quantize_midrise(u, bits, y_min, y_max)
    d = y_max - y_min
    y1 = clamp(u, y_min, y_max)
    y2 = (y1 - y_min) / d # between 0 and 1
    Δ = 2^Int(bits)-1
    y3 = round(y2 * Δ) / Δ # quantized between 0 and 1
    y4 = y3*d + y_min
    return y4
end

function quantize_midtread(u, bits, y_min, y_max)
    Δ = (y_max - y_min) / (2^Int(bits)-1)
    # clamp(Δ * floor(u / Δ + 0.5), y_min, y_max)
    k = sign(u) * max(0, floor((abs(u) - Δ/2) / Δ + 1))
    y0 = sign(k) * (Δ/2 + Δ*(abs(k)-1/2))
    y1 = iszero(y0) ? zero(y0) : y0 # remove -0.0
    return clamp(y1, y_min, y_max - Δ/2)
end

function quantize(u, bits, y_min, y_max, midrise)
    midrise ? quantize_midrise(u, bits, y_min, y_max) : quantize_midtread(u, bits, y_min, y_max)
end

@register_symbolic quantize(u::Real, bits::Real, y_min::Real, y_max::Real, midrise::Bool)

"""
    ExponentialFilter(a = 0.1)

Exponential filtering (first-order filter) with input-output relation ``y(z) = (1 - a) y(z-1) + a u(z-1)``, transfer function
```math
Y(z) = \\dfrac{a}{1 - (1 - a) z^{-1}} U(z)
```

# Parameters:
- `a`: Filter parameter `[0, 1]`, a small value implies stronger filtering. 

# Variables:
- `u`: Input signal
- `y`: Output signal

# Connectors:
- `input::RealInput`: Input signal
- `output::RealOutput`: Output signal
"""
@mtkmodel ExponentialFilter begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        z = ShiftIndex()
    end
    @parameters begin
        a = 0.1, [description = "Filter parameter"]
    end
    @equations begin
        y(z) ~ (1 - a) * y(z-1) + a * u(z)
    end
end

"""
    MovingAverageFilter(N = 3)

Exponential filtering with input-output relation ``y(z) = sum(u(z-i) for i in 0:N-1) / N``.

Please note: this implementation of a moving average filter is not optimized for very large number of filter taps `N`.

# Parameters:
- `N`: (structural) Number of samples to average over
# Variables:
- `u`: Input signal
- `y`: Output signal

# Connectors:
- `input::RealInput`: Input signal
- `output::RealOutput`: Output signal
"""
@mtkmodel MovingAverageFilter begin
    @extend u, y = siso = SISO()
    @structural_parameters begin
        z = ShiftIndex()
        N = 3
    end
    @equations begin
        y(z) ~ sum(u(z-i) for i in 0:N-1) / N
    end
end

"""
    DiscreteOnOffController(b = 0.1, bool = true)

Discrete-time On-Off controller with hysteresis. The controller switches between two states based on the error signal `reference-input`. The controller is in the on-state if the error signal is within the bandwidth `b` around the reference signal, and in the off-state otherwise.

# Connectors:
- `reference`: The reference signal to the controller
- `input`: The measurement feedback
- `output`: The control signal output

# Parameters:
- `b`: Bandwidth around reference signal within which the controller does not react
- `bool`: (structural) If true (default), the controller switches between 0 and 1. If false, the controller switches between -1 and 1.
- `k`: Controller gain. The output of the controller is scaled by this gain, i.e., `k = 2, bool = false` will result in an output of -2 or 2.
"""
@mtkmodel DiscreteOnOffController begin
    @extend u, y = siso = SISO()
    @components begin
        reference = RealInput()
    end
    @structural_parameters begin
        z = ShiftIndex()
        bool = true
    end
    @parameters begin
        b = 0.1, [description = "Bandwidth around reference signal"]
        k = 1, [description = "Controller gain"]
    end
    @variables begin
        s(t)=true, [description = "Internal variable"]
    end
    @equations begin
        s(z) ~ (y(z-1) == k) & (u(z) < reference.u(z) + b/2) | (u(z) < reference.u(z) - b/2)
        if bool
            y(z) ~ k*s(z)
        else
            y(z) ~ k*(2*s(z) - 1)
        end
    end
end

"""
    SampleWithADEffects(quantized = false, noisy = false)

A sampler with additional effects that appear in practical systems, such as measurement noise and quantization.

The operations occur in the order
1. Sampling
2. Noise addition
3. Quantization

# Structural parameters:
- `quantized`: If true, the output is quantized. When this option is used, the output is quantized to the number of bits specified by the `bits` parameter. The quantization is midrise if `midrise = true`, otherwise it is midtread. The output is also limited to the range `[y_min, y_max]`.
- `noisy`: If true, the output is corrupted by additive white Gaussian noise with standard deviation `sigma` (defaults to 0.1). If `noisy = false`, the noise block is a unit gain.
- `dt`: Sample interval of the sampler. If not specified, the sample interval is inferred from the clock of the system.
- `clock`: Clock signal of the system. If not specified, the sample interval is inferred from the clock of the system. If `clock` is specified, the parameter `dt` has no effect.

# Parameters:
- `y_min`: Lower limit of output, defaults to -1. Only used if `quantized = true`.
- `y_max`: Upper limit of output, defaults to 1. Only used if `quantized = true`.
- `bits`: Number of bits of quantization, defaults to 8 (256 output levels between `y_min` and `y_max`). Only used if `quantized = true`.
- `sigma`: Standard deviation of the additive Gaussian noise, defaults to 0.1. Only used if `noisy = true`.
"""
@mtkmodel SampleWithADEffects begin
    @extend input, output = siso = SISO()
    @structural_parameters begin
        dt = nothing
        clock = (dt === nothing ? InferredDiscrete() : Clock(dt))
        midrise = true
        quantized = false
        noisy = false
    end
    @parameters begin
        y_min = -1, [description = "Lower limit of output"]
        y_max = 1, [description = "Upper limit of output"]
        bits = 8, [description = "Number of bits of quantization"]
        sigma = 0.1, [description = "Standard deviation of the additive noise"]
    end
    @components begin
        sampler = Sampler(; clock)
        if noisy
            noise = NormalNoise(; sigma, additive = true)
        else
            noise = Gain(; k = 1)
        end
        quantization = Quantization(; bits, y_min, y_max, midrise, quantized)
    end
    @equations begin
        connect(input, sampler.input)
        connect(sampler.output, noise.input)
        connect(noise.output, quantization.input)
        connect(quantization.output, output)
    end
end
