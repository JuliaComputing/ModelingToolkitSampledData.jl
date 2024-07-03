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
            kp = 2, ki = 2, Imethod = :forward, with_D = false)
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
model = structural_simplify(IRSystem(model))

Tf = 5
timevec = 0:(dt):Tf

import ControlSystemsBase as CS
import ControlSystemsBase: c2d, tf, feedback, lsim
P = CS.c2d(CS.ss([-1], [1], [1], 0), dt)
C = CS.c2d(CS.ss([0], [1], [2], [2]), dt, :fwdeuler)

# Test the output of the continuous partition
G = feedback(P * C)
res = lsim(G, (x, t) -> [0.5], timevec)
y = res.y[:]

prob = ODEProblem(model,
    [model.plant.x => 0.0; model.controller.kp => 2.0; model.controller.ki => 2.0;
     model.controller.eI => 0.0],
    (0.0, Tf))

sol = solve(prob,
    Tsit5(),
    abstol = 1e-8,
    reltol = 1e-8)

# plot(timevec, [y sol(timevec, idxs = model.plant.output.u)[:]], m = :o, lab = ["CS" "MTK"])
# display(current())

@test sol(timevec, idxs = model.plant.output.u)[:]≈y rtol=1e-6

@test_skip begin
    # Test the output of the discrete partition
    G = feedback(C, P)
    res = lsim(G, (x, t) -> [0.5], timevec)
    y = res.y[:]
    @test_broken sol(timevec .+ 1e-10, idxs = model.controller.output.u)≈y rtol=1e-8 # Broken due to discrete observed
    # plot([y sol(timevec .+ 1e-12, idxs=model.controller.output.u)], lab=["CS" "MTK"])
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
            K = 2, Ti = 1, Imethod = :forward, with_D = false)
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
model = structural_simplify(IRSystem(model))

Tf = 5
timevec = 0:(dt):Tf

import ControlSystemsBase as CS
import ControlSystemsBase: c2d, tf, feedback, lsim
P = CS.c2d(CS.ss([-1], [1], [1], 0), dt)
C = CS.c2d(CS.ss([0], [1], [2], [2]), dt, :fwdeuler)

# Test the output of the continuous partition
G = feedback(P * C)
res = lsim(G, (x, t) -> [0.5], timevec)
y = res.y[:]

prob = ODEProblem(model,
    [model.plant.x => 0.0; model.controller.eI => 0.0],
    (0.0, Tf))

sol = solve(prob,
    Tsit5(),
    abstol = 1e-8,
    reltol = 1e-8)

# plot(timevec, [y sol(timevec, idxs = model.plant.output.u)[:]], m = :o, lab = ["CS" "MTK"])
# display(current())

@test sol(timevec, idxs = model.plant.output.u)[:]≈y rtol=1e-6

@test_skip begin
    # Test the output of the discrete partition
    G = feedback(C, P)
    res = lsim(G, (x, t) -> [0.5], timevec)
    y = res.y[:]
    @test_broken sol(timevec .+ 1e-10, idxs = model.controller.output.u)≈y rtol=1e-8 # Broken due to discrete observed
    # plot([y sol(timevec .+ 1e-12, idxs=model.controller.output.u)], lab=["CS" "MTK"])
end

# ==============================================================================
## Delay
# ==============================================================================

@testset "delay" begin
    @mtkmodel DelayModel begin
        @components begin
            fake_plant = FirstOrder(T = 1e-4) # Included due to bug with only discrete-time systems
            input = Step(start_time = 2, smooth = false)
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
            input = Step(start_time = 2, smooth = false)
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

@testset "NormalNoise" begin
    k = ShiftIndex(Clock(1))

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
    prob = ODEProblem(ssys, [], (0.0, 10.0))
    sol = solve(prob, Tsit5())
    # TODO: add tests
end
