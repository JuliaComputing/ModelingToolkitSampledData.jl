using ModelingToolkitSampledData
using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq: ReturnCode.Success
using JuliaSimCompiler
using Test
Difference = ModelingToolkitSampledData.Difference

#=
Testing strategy:
The general strategy is to test systems using simple inputs where the solution
is known on closed form. For algebraic systems (without differential variables),
an integrator with a constant input is often used together with the system under test.
=#

@testset "only discrete" begin
    @info "Testing only discrete"
    dt = 0.5
    c = Clock(dt)
    k = ShiftIndex(c)

    @mtkmodel TestDiscreteOnly begin
        @variables begin
            x(t) = 1
        end
        @equations begin
            x(k) ~ SampleTime()*x(k-1)
        end
    end

    @named model = TestDiscreteOnly()
    model = complete(model)
    ssys = structural_simplify(IRSystem(model))
    prob = ODEProblem(ssys, [model.x(k-1) => 1.0], (0.0, 10.0))
    sol = solve(prob, Tsit5())
    @test sol[model.x] == dt .^ (1:21)
end

@testset "Integrator" begin
    clock = Clock(0.1)
    k = ShiftIndex(clock)
    @named c = Constant(; k = 1)
    @named sampler = Sampler(; clock)
    @named intc = Integrator()

    # Backward
    @named int = DiscreteIntegrator(x = 1)
    @named iosys = ODESystem(
        [connect(c.output, sampler.input)
         connect(sampler.output, int.input)
         connect(c.output, intc.input)],
        t,
        systems = [sampler, int, intc, c])
    model = complete(iosys)
    sys = structural_simplify(IRSystem(model))
    prob = ODEProblem(sys, Pair[int.x(k - 1) => 1
                                int.u(k - 1) => 0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test_skip sol.prob.kwargs[:disc_saved_values][1].t ≈ 0:sampletime(clock):1
    @test sol[model.int.x] ≈ range(1.1, step = sampletime(clock), length = 11)

    # Forward
    @named int = DiscreteIntegrator(x = 1, method = :forward)
    @named iosys = ODESystem(
        [connect(c.output, sampler.input)
         connect(sampler.output, int.input)
         connect(c.output, intc.input)],
        t,
        systems = [sampler, int, intc, c])
    model = complete(iosys)
    sys = structural_simplify(IRSystem(model))
    prob = ODEProblem(sys, Pair[int.x(k - 1) => 1
                                int.u(k - 1) => 0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test_skip sol.prob.kwargs[:disc_saved_values][1].t ≈ 0:sampletime(clock):1
    @test sol[model.int.x] ≈
          range(1.0, step = sampletime(clock), length = 11)

    # Tustin
    @named int = DiscreteIntegrator(x = 1, method = :tustin)
    @named iosys = ODESystem(
        [connect(c.output, sampler.input)
         connect(sampler.output, int.input)
         connect(c.output, intc.input)],
        t,
        systems = [sampler, int, intc, c])
    model = complete(iosys)
    sys = structural_simplify(IRSystem(model))
    prob = ODEProblem(sys, Pair[int.x(k - 1) => 1
                                int.u(k - 1) => 0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test_skip sol.prob.kwargs[:disc_saved_values][1].t ≈ 0:sampletime(clock):1
    @test sol[model.int.x] ≈
          range(1.05, step = sampletime(clock), length = 11)
end


# ==============================================================================
## DiscretePIDParallel
# ==============================================================================

dt = 0.05
k = ShiftIndex()

@mtkmodel ClosedLoop begin
    @components begin
        plant = FirstOrder(k = 1, T = 1)
        sampler = Sampler(; dt)
        zoh = ZeroOrderHold()
        controller = DiscretePIDParallel(
            kp = 2, ki = 2, Imethod = :backward, with_D = false) # NOTE: not sure why tests pass with backward euler here but fwdeuler in ControlSystemsBase
        ref = Constant(k = 0.5)
    end
    @equations begin
        connect(ref.output, controller.reference)
        connect(controller.ctr_output, zoh.input)
        connect(zoh.output, plant.input)
        connect(plant.output, sampler.input)
        connect(sampler.output, controller.measurement)
    end
end

#
@named model = ClosedLoop()
model = complete(model)
# ci, varmap = infer_clocks(expand_connections(model))
ssys = structural_simplify(IRSystem(model))

Tf = 5
timevec = 0:(dt):Tf

import ControlSystemsBase as CS
let (; c2d, tf, feedback, lsim) = CS
    
    P = CS.c2d(CS.ss([-1], [1], [1], 0), dt)
    # C = CS.c2d(CS.ss([0], [1], [2], [2]), dt, :fwdeuler)
    C = CS.ss([1], [1], [2*dt], [2], dt)

    # Test the output of the continuous partition
    G = feedback(P * C)
    res = lsim(G, (x, t) -> [0.5], timevec)
    y = res.y[:]

    prob = ODEProblem(ssys,
        [model.plant.x => 0.0; model.controller.kp => 2.0; model.controller.ki => 2.0;
        model.controller.eI(k-1) => 0.0; model.controller.I(k-1) => 0.0],
        (0.0, Tf))

    sol = solve(prob,
        Tsit5(),
        abstol = 1e-8,
        reltol = 1e-8)

    # plot(timevec, [y sol(timevec, idxs = model.plant.output.u)[:]], m = :o, lab = ["CS" "MTK"])
    # display(current())

    @test sol(timevec, idxs = model.plant.output.u)[:]≈y rtol=1e-5

    # Test the output of the discrete partition
    G = feedback(C, P)
    res = lsim(G, (x, t) -> [0.5], timevec)
    y = res.y[:]
    @test_broken sol(timevec, idxs = model.plant.input.u)≈y rtol=1e-8 
    # plot([y sol(timevec, idxs=model.plant.input.u).u], lab=["CS" "MTK"])
    # plot(timevec, [y sol(timevec, idxs = model.controller.u)[:]], m = :o, lab = ["CS" "MTK"])
end
# ==============================================================================
## DiscretePIDStandard
# ==============================================================================

dt = 0.05
k = ShiftIndex()

@mtkmodel ClosedLoop begin
    @components begin
        plant = FirstOrder(k = 1, T = 1)
        sampler = Sampler(; dt)
        zoh = ZeroOrderHold()
        controller = DiscretePIDStandard(
            K = 2, Ti = 1, Imethod = :backward, with_D = false)
        ref = Constant(k = 0.5)
    end
    @equations begin
        connect(ref.output, controller.reference)
        connect(controller.ctr_output, zoh.input)
        connect(zoh.output, plant.input)
        connect(plant.output, sampler.input)
        connect(sampler.output, controller.measurement)
    end
end

#
@named model = ClosedLoop()
model = complete(model)
# ci, varmap = infer_clocks(expand_connections(model))
ssys = structural_simplify(IRSystem(model))

Tf = 5
timevec = 0:(dt):Tf

import ControlSystemsBase as CS
let (; c2d, tf, feedback, lsim) = CS
    P = CS.c2d(CS.ss([-1], [1], [1], 0), dt)
    C = CS.c2d(CS.ss([0], [1], [2], [2]), dt, :fwdeuler)

    # Test the output of the continuous partition
    G = feedback(P * C)
    res = lsim(G, (x, t) -> [0.5], timevec)
    y = res.y[:]

    prob = ODEProblem(ssys,
        [model.plant.x => 0.0; model.controller.eI(k-1) => 0.0; model.controller.I(k-1) => 0.0],
        (0.0, Tf))

    sol = solve(prob,
        Tsit5(),
        abstol = 1e-8,
        reltol = 1e-8)

    # plot(timevec, [y sol(timevec, idxs = model.plant.output.u)[:]], m = :o, lab = ["CS" "MTK"])
    # display(current())

    @test sol(timevec, idxs = model.plant.output.u)[:]≈y rtol=1e-6

    # Test the output of the discrete partition
    G = feedback(C, P)
    res = lsim(G, (x, t) -> [0.5], timevec)
    y = res.y[:]
    @test_broken sol(timevec, idxs = model.controller.u)≈y rtol=1e-8 # Broken due to LLVM error: ERROR: KeyError: key shift(controller₊I, -1, false, nothing) not found
    # plot([y sol(timevec .+ 1e-12, idxs=model.controller.u).u], lab=["CS" "MTK"])
end
# ==============================================================================
## Delay
# ==============================================================================

@testset "delay" begin
    @mtkmodel DelayModel begin
        @components begin
            fake_plant = FirstOrder(T = 1e-4) # Included due to bug with only discrete-time systems
            input = Step(start_time = 2-1e-9, smooth = false) # shift by 1e-9 to account for MTKStdlib using > instead of >=
            sampler = Sampler(; dt = 1)
            delay = Delay(n = 3)
            zoh = ZeroOrderHold()
        end
        @equations begin
            connect(input.output, sampler.input)
            connect(sampler.output, delay.input)
            connect(delay.output, zoh.input)
            connect(zoh.output, fake_plant.input)
        end
    end

    @named m = DelayModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(
        ssys, [m.delay.u(k - 3) => 0, m.delay.u(k - 2) => 0, m.delay.u(k - 1) => 0], (
            0.0, 10.0))
    sol = solve(prob, Tsit5())

    @test reduce(vcat, sol((0:10) .+ 1e-2))[:]≈[zeros(5); ones(6)] atol=1e-2
end

# ==============================================================================
## Difference
# ==============================================================================
using ModelingToolkitStandardLibrary.Blocks

@testset "Difference" begin
    k = ShiftIndex(Clock(1))

    @mtkmodel DiffModel begin
        @components begin
            input = Step(start_time = 2-1e-9, smooth = false)
            diff = Difference(z = k)
            zoh = ZeroOrderHold()
            plant = FirstOrder(T = 1e-4) # Included due to bug with only discrete-time systems
        end
        @equations begin
            connect(input.output, diff.input)
            connect(diff.output, zoh.input)
            connect(zoh.output, plant.input)
        end
    end

    @named m = DiffModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, Dict(m.diff.u(k - 1) => 0), (0.0, 10.0))
    sol = solve(prob, Tsit5(), dtmax = 0.01)
    @test reduce(vcat, sol((0:10) .+ 1e-2))[:]≈[zeros(2); 1; zeros(8)] atol=1e-2
end


# ==============================================================================
## Noise sources
# ==============================================================================
using ModelingToolkitStandardLibrary.Blocks
using Statistics

@testset "NormalNoise" begin
    k = ShiftIndex(Clock(0.01))

    @mtkmodel NoiseModel begin
        @components begin
            noise = NormalNoise(z = k)
            zoh = ZeroOrderHold(z = k)
            plant = FirstOrder(T = 1e-4) # Included due to bug with only discrete-time systems
        end
        @equations begin
            connect(noise.output, zoh.input)
            connect(zoh.output, plant.input)
        end
    end

    @named m = NoiseModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, [m.noise.y(k-1) => 0], (0.0, 10.0))
    sol = solve(prob, Tsit5())
    @test !all(iszero, sol.u)
    tv = 0:k.clock.dt:sol.t[end]
    @test std(sol(tv, idxs = m.noise.y)) ≈ 1 rtol=0.1
    @test mean(sol(tv, idxs = m.noise.y)) ≈ 0 atol=0.08
end

@testset "UniformNoise" begin
    k = ShiftIndex(Clock(0.01))

    @mtkmodel NoiseModel begin
        @components begin
            noise = UniformNoise(z = k)
            zoh = ZeroOrderHold(z = k)
            plant = FirstOrder(T = 1e-4) # Included due to bug with only discrete-time systems
        end
        @equations begin
            connect(noise.output, zoh.input)
            connect(zoh.output, plant.input)
        end
    end

    @named m = NoiseModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, [m.noise.y(k-1) => 0], (0.0, 10.0))
    sol = solve(prob, Tsit5())
    @test !all(iszero, sol.u)
    tv = 0:k.clock.dt:sol.t[end]
    @test minimum(sol(tv, idxs = m.noise.y)) ≈ 0 atol=0.02
    @test maximum(sol(tv, idxs = m.noise.y)) ≈ 1 atol=0.02
end


# https://github.com/SciML/ModelingToolkit.jl/issues/2843
# @testset "StateSpace" begin
#     k = ShiftIndex(Clock(1))

#     @mtkmodel PlantModel begin
#         @components begin
#             input = Constant(k=1)
#             plant = DiscreteStateSpace(z = k, A=[1;;], B=[1;;], C=[1;;], D=[0;;])
#         end
#         @variables begin
#             x(t) # Dummy variable
#         end
#         @equations begin
#             connect(input.output, plant.input)
#             D(x) = 1
#         end
#     end

#     @named m = PlantModel()
#     m = complete(m)
#     ssys = structural_simplify(IRSystem(m))
#     prob = ODEProblem(ssys, Dict(m.plant.u(k - 1) => 0), (0.0, 10.0))
#     sol = solve(prob, Tsit5(), dtmax = 0.01)
#     @test reduce(vcat, sol((0:10) .+ 1e-2))[:]≈[zeros(2); 1; zeros(8)] atol=1e-2
# end*

##
@testset "SlewRateLimiter" begin
    @info "Testing SlewRateLimiter"
    cl = Clock(0.1)
    z = ShiftIndex(cl)
    @mtkmodel SlewRateLimiterModel begin
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
    @named m = SlewRateLimiterModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, [m.limiter.y(z-1) => 0], (0.0, 2.0))
    sol = solve(prob, Tsit5(), dtmax=0.01)
    # plot(sol, idxs=[m.input.output.u, m.limiter.y], title="Slew rate limited sine wave")
    @test maximum(diff(sol[m.limiter.y])) ≈ 0.4
    @test minimum(diff(sol[m.limiter.y])) ≈ -0.3
end

##
@testset "quantization" begin
    @info "Testing quantization"

    function test_quant(y_min, y_max, bits)
        u = y_min:(1/2^bits):y_max
        y = ModelingToolkitSampledData.quantize_midrise.(u, bits, y_min, y_max)
        uy = unique(y)
        @test uy ≈ range(y_min, stop=y_max, length=2^bits)
    end

    test_quant(-1, 1, 2) # Symmetric
    test_quant(-1, 2, 2) # Not symmetric
    test_quant(-1, 1, 3) # Symmetric, uneven number of bits
    test_quant(-5, -2, 2) # Only negative
    test_quant(5, 12, 2) # Only positive
    test_quant(-5, 12, 20) # Large number of bits



    function test_quant2(y_min, y_max, bits)
        u = y_min:(1/2^bits):y_max
        y = ModelingToolkitSampledData.quantize_midtread.(u, bits, y_min, y_max)
        uy = unique(y) 
        # @test (2^bits - 2 <= length(uy) <= 2^bits) # This check is not reliable since there might be one ulp difference when the last clamp is applied
        @test maximum(y) <= y_max
        @test minimum(y) >= y_min
    end

    test_quant2(-1, 1, 2) # Symmetric
    test_quant2(-1, 2, 2) # Not symmetric
    test_quant2(-1, 1, 3) # Symmetric, uneven number of bits
    test_quant2(-5, -2, 2) # Only negative
    test_quant2(5, 12, 2) # Only positive
    test_quant2(-5, 12, 20) # Large number of bits

    z = ShiftIndex(Clock(0.1))
    @mtkmodel QuantizationModel begin
        @components begin
            input = Sine(amplitude=1, frequency=1)
            quant = Quantization(; z, bits=2)
        end
        @equations begin
            connect(input.output, quant.input)
        end
    end
    @named m = QuantizationModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, [], (0.0, 10.0))
    sol = solve(prob, Tsit5(), dtmax=0.01)
    y = sol[m.quant.y]
    uy = unique(y)
    @test length(uy) == 4


    @mtkmodel QuantizationModel2 begin
        @components begin
            input = Sine(amplitude=1, frequency=1)
            quant = Quantization(; z, bits=2, midrise=false)
        end
        @equations begin
            connect(input.output, quant.input)
        end
    end
    @named m = QuantizationModel2()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, [], (0.0, 10.0))
    sol = solve(prob, Tsit5(), dtmax=0.01)
    y = sol[m.quant.y]
    uy = unique(y)
    @test length(uy) <= 4
    @test 0 ∈ uy
end

@testset "OnOff" begin
    @info "Testing OnOff"
    cl = Clock(0.1)
    z = ShiftIndex(cl)
    @mtkmodel OnOffModel begin
        @components begin
            onoff = DiscreteOnOffController(; z, bool=false)
            c = Constant(k=1)
        end
        @variables begin
            x(t) = 0 
        end
        @equations begin
            onoff.u ~ Sample(cl)(x)
            connect(c.output, onoff.reference)
            D(x) ~ 0.1x + Hold(onoff.y)
        end
    end
    @named m = OnOffModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, [m.onoff.y(z-1) => 0], (0.0, 4.0))
    sol = solve(prob, Tsit5(), dtmax=0.1)
    # plot(sol, idxs=[m.x, m.onoff.y], title="On-off control of an unstable first-order system")
    @test 0.89 <= sol(4, idxs=m.x) <= 1.11
end

@testset "ExponentialFilter" begin
    @info "Testing ExponentialFilter"
    z = ShiftIndex(Clock(0.1))
    @mtkmodel ExponentialFilterModel begin
        @components begin
            input = Step(start_time=1, smooth=false)
            filter = ExponentialFilter(; a = 0.1, z)
        end
        @variables begin
            x(t) = 0 # Dummy variable to workaround JSCompiler bug
        end
        @equations begin
            connect(input.output, filter.input)
            D(x) ~ 0
        end
    end
    @named m = ExponentialFilterModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, [m.filter.y(z-1) => 0], (0.0, 10.0))
    sol = solve(prob, Tsit5(), dtmax=0.1)
    @test sol(10, idxs=m.filter.y) ≈ 1 atol=0.001
    @test sol(0.999, idxs=m.filter.y) == 0
    @test sol(1.1, idxs=m.filter.y) > 0

end

@testset "MovingAverageFilter" begin
    @info "Testing MovingAverageFilter"
    z = ShiftIndex(Clock(0.1))
    @mtkmodel MovingAverageFilterModel begin
        @components begin
            input = Step(start_time=1, smooth=false)
            filter = MovingAverageFilter(; N=3, z)
        end
        @variables begin
            x(t) = 0 # Dummy variable to workaround JSCompiler bug
        end
        @equations begin
            connect(input.output, filter.input)
            D(x) ~ 0
        end
    end

    @named m = MovingAverageFilterModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, [m.filter.u(z-i) => 0 for i = 0:3], (0.0, 2.0))
    sol = solve(prob, Tsit5(), dtmax=0.1)
    # plot(sol, idxs=m.filter.y)
    @test sol(1.5, idxs=m.filter.y) == 1
    @test sol(0.999, idxs=m.filter.y) == 0
    @test_broken 0 < sol(1.1, idxs=m.filter.y) < 1

    @test_broken count(!∈((1,0)), sol[m.filter.y]) == 2 # With 3 taps there are two steps where tap values have both 0 and 1
end

@testset "sampling with AD effects" begin
    @info "Testing sampling with AD effects"
    z = ShiftIndex()
    @mtkmodel PracticalSamplerModel begin
        @components begin
            input = Sine(amplitude=1.2, frequency=1, smooth=false)
            sampling = SampleWithADEffects(; dt=0.03, bits=3, y_min = -1, y_max = 1, sigma = 0.1, noisy = true, quantized=true, midrise=true)
        end
        @variables begin
            x(t) = 0 # Dummy variable to work around a bug for models without continuous-time state
        end
        @equations begin
            connect(input.output, sampling.input)
            D(x) ~ Hold(sampling.y) # Dummy equation
        end
    end
    @named m = PracticalSamplerModel()
    m = complete(m)
    ssys = structural_simplify(IRSystem(m))
    prob = ODEProblem(ssys, [m.sampling.noise.y(z-1) => 0], (0.0, 2.0))
    sol = solve(prob, Tsit5())
    
    @test length(unique(sol[m.sampling.y])) == 8
    @test maximum(abs, sol(0:0.03:1, idxs=m.sampling.y) - sol(0:0.03:1, idxs=m.sampling.u)) < 0.3

end
