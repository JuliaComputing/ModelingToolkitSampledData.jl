# Measurement noise and corruption
Measurement noise is practically always present in signals originating from real-world sensors. In a sampled-data system, analyzing the influence of measurement noise using simulation is relatively straight forward. Below, we add Gaussian white noise to the speed sensor signal in the DC motor example. The noise is added using the [`NormalNoise`](@ref) block.

This block has two modes of operation
1. If `additive = false` (default), the block has the connector `output` only, and this output is the noise signal.
2. If `additive = true`, the block has the connectors `input` and `output`, and the output is the sum of the input and the noise signal, i.e., the noise is _added_ to the input signal. This mode makes it convenient to add noise to a signal in a sampled-data system.

## Example: Noise
```@example NOISE
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitSampledData
using JuliaSimCompiler
using OrdinaryDiffEq
using Plots

R = 0.5 # [Ohm] armature resistance
L = 4.5e-3 # [H] armature inductance
k = 0.5 # [N.m/A] motor constant
J = 0.02 # [kg.m²] inertia
f = 0.01 # [N.m.s/rad] friction factor
tau_L_step = -0.3 # [N.m] amplitude of the load torque step
nothing # hide

z = ShiftIndex()

@mtkmodel NoisyClosedLoop begin
    @components begin
        ground = Ground()
        source = Voltage()
        ref = Blocks.Step(height = 1, start_time = 0, smooth = false)
        sampler = Sampler(dt = 0.002)
        noise = NormalNoise(sigma = 0.1, additive = true)
        pi_controller = DiscretePIDStandard(
            K = 1, Ti = 0.035, u_max = 10, with_D = false)
        zoh = ZeroOrderHold()
        R1 = Resistor(R = R)
        L1 = Inductor(L = L)
        emf = EMF(k = k)
        fixed = Fixed()
        load = Torque()
        load_step = Blocks.Step(height = tau_L_step, start_time = 1.3)
        inertia = Inertia(J = J)
        friction = Damper(d = f)
        speed_sensor = SpeedSensor()
        angle_sensor = AngleSensor()
    end

    @equations begin
        connect(fixed.flange, emf.support, friction.flange_b)
        connect(emf.flange, friction.flange_a, inertia.flange_a)
        connect(inertia.flange_b, load.flange)
        connect(inertia.flange_b, speed_sensor.flange, angle_sensor.flange)
        connect(load_step.output, load.tau)
        connect(ref.output, pi_controller.reference)
        connect(speed_sensor.w, sampler.input)
        connect(sampler.output, noise.input)
        connect(noise.output, pi_controller.measurement)
        connect(pi_controller.ctr_output, zoh.input)
        connect(zoh.output, source.V)
        connect(source.p, R1.p)
        connect(R1.n, L1.p)
        connect(L1.n, emf.p)
        connect(emf.n, source.n, ground.g)
    end
end


@named noisy_model = NoisyClosedLoop()
noisy_model = complete(noisy_model)
ssys = structural_simplify(IRSystem(noisy_model)) # Conversion to an IRSystem from JuliaSimCompiler is required for sampled-data systems

noise_prob = ODEProblem(ssys, [unknowns(noisy_model) .=> 0.0; noisy_model.pi_controller.I(z-1) => 0; noisy_model.pi_controller.eI(z-1) => 0; noisy_model.noise.y(z-1) => 0], (0, 2.0))
noise_sol = solve(noise_prob, Tsit5())

figy = plot(noise_sol, idxs=noisy_model.noise.y, label = "Measured speed", )
plot!(noise_sol, idxs=noisy_model.inertia.w, ylabel = "Angular Vel. [rad/s]",
    label = "Actual speed", legend=:bottomleft, dpi=600, l=(2, :blue))
figu = plot(noise_sol, idxs=noisy_model.source.V.u, label = "Control signal [V]", )
plot(figy, figu, plot_title = "DC Motor with Discrete-time Speed Controller")
```

## Linear analysis of noise
Propagation of Gaussian noise through linear time-invariant systems is well understood, the stationary covariance of the output can be computed by solving a Lyapunov equation. Unfortunately, ModelingToolkit models that contain both continuous time and discrete time components cannot yet be linearized and linear analysis is thus made slightly harder. Below, we instead use a data-driven linearization approach where we use recorded signals from the simulation and fit a linear model using subspace-based identification. The function `subspaceid` below is provided by the package [ControlSystemIdentification.jl](https://baggepinnen.github.io/ControlSystemIdentification.jl/stable/).

We let the angular velocity of the inertia be the output, and the output of the noise block as well as the output of the load disturbance be the inputs. 

```@example NOISE
using ControlSystemIdentification, ControlSystemsBase
Tf = 20
prob2 = remake(noise_prob, p=Dict(noisy_model.load_step.height=>0.0), tspan=(0.0, Tf))
noise_sol = solve(prob2, Tsit5())
tv = 0:0.002:Tf
y = noise_sol(tv, idxs=noisy_model.inertia.w) |> vec
un = noise_sol(tv, idxs=noisy_model.noise.y)-y |> vec
ud = noise_sol(tv, idxs=noisy_model.load_step.output.u) |> vec
d = iddata(y', [un ud]', 0.002)
lsys,_ = newpem(d, 4, focus=:simulation, zeroD=false)
```
With an LTI model available, we can ask for the theoretical output covariance we should obtain if we feed a white noise signal with covariance matrix ``0.1^2 I`` through the noise input of the system. We compare this to the actual output covariance obtained from the simulation (discarding the initial transient as well as the transient caused by the load disturbance).
```@example NOISE
sqrt(covar(lsys[1,1],0.1^2*I)), std(y[[50:648; 750:end]])
```

## Noise filtering
No discrete-time filter components are available yet. You may, e.g.
- Add exponential filtering using `xf(k) ~ (1-α)xf(k-1) + α*x(k)`, where `α` is the filter coefficient and `x` is the signal to be filtered.
- Add moving average filtering using `xf(k) ~ 1/N sum(i->x(k-i), i=0:N-1)`, where `N` is the number of samples to average over.

## Colored noise
Colored noise can be achieved by filtering white noise through a filter with the desired spectrum. No components are available for this yet.

## Internal details
Internally, a random number generator from [StableRNGs.jl](https://github.com/JuliaRandom/StableRNGs.jl) is used to produce reproducible streams of random numbers. Each draw of a random number is seeded by `hash(t, hash(seed))`, where `seed` is a parameter in the noise source component, and `t` is the current simulation time. This ensures that
1. The user can alter the stream of random numbers with `seed`.
2. Multiple calls to the random number generator at the same time step all return the same number.

## Quantization

A signal may be quantized to a fixed number of levels (e.g., 8-bit) using the [`Quantization`](@ref) block. This may be used to simulate, e.g., the quantization that occurs in a AD converter. Below, we have a simple example where a sine wave is quantized to 2 bits (4 levels), limited between -1 and 1:
```@example QUANT
using ModelingToolkit, ModelingToolkitSampledData, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D
z = ShiftIndex(Clock(0.1))
@mtkmodel QuantizationModel begin
    @components begin
        input = Sine(amplitude=1.5, frequency=1)
        quant = Quantization(; z, bits=2, y_min = -1, y_max = 1)
    end
    @variables begin
        x(t) = 0 # Dummy variable to work around a bug for models without continuous-time state
    end
    @equations begin
        connect(input.output, quant.input)
        D(x) ~ 0 # Dummy equation
    end
end
@named m = QuantizationModel()
m = complete(m)
ssys = structural_simplify(IRSystem(m))
prob = ODEProblem(ssys, [], (0.0, 2.0))
sol = solve(prob, Tsit5())
plot(sol, idxs=m.input.output.u)
plot!(sol, idxs=m.quant.y, label="Quantized output")
```



### Different quantization modes
With the default option `midrise = true`, the output of the quantizer is always between `y_min` and `y_max` inclusive, and the number of distinct levels it can take is `2^bits`. The possible values are given by
```@example
bits = 2; y_min = -1; y_max = 1
collect(range(y_min, stop=y_max, length=2^bits))
```
Notably, these possible levels _do not include 0_. If `midrise = false`, a mid-tread quantizer is used instead. The two options are visualized below:
```@example QUANT
y_min = -1; y_max = 1; bits = 2
u = y_min:0.01:y_max
y_mr = ModelingToolkitSampledData.quantize_midrise.(u, bits, y_min, y_max)
y_mt = ModelingToolkitSampledData.quantize_midtread.(u, bits, y_min, y_max)
plot(u, [y_mr y_mt], label=["Midrise" "Midtread"], xlabel="Input", ylabel="Output", framestyle=:zerolines, l=2, seriestype=:step)
```
Note how the default mid-rise quantizer mode has a rise at the middle of the interval, while the mid-tread mode has a flat region (a tread) centered around the middle of the interval.

The default option `midrise = true` includes both end points as possible output values, while `midrise = false` does not include the upper limit.
