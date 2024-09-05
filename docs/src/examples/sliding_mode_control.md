# Sliding-mode control
This example demonstrates how to model a system with a discrete-time sliding-mode controller using ModelingToolkit. The system consists of a second-order plant with a disturbance and a super-twisting sliding-mode controller. The controller is implemented as a discrete-time system, and the plant is modeled as a continuous-time system. 

When designing an SMC controller, we must choose a _switching function_ ``σ`` that produces the sliding variable ``s = σ(x, p, t)``. The sliding surface ``s=0`` must be chosen such that the sliding variable ``s`` exhibits desireable properties, i.e., converges to the desired state with stable dynamics. The control law is chosen to drive the system from any state ``x : s ≠ 0`` to the sliding surface ``s=0``. The sliding surface is commonly chosen as an asymptotically stable system with order equal to the ``n_x-n_u``, where ``n_x`` is the dimension of the state in the system to be controlled, and ``n_u`` is the number of inputs.

Since the dynamics in this example is of relative degree ``r=2``, we will choose a switching surface corresponding to a stable first-order system (``r-1=1``). We will choose the system
```math
ė = -e
```
which yields the switching variable ``s = ė + e``, encoded in the function ``s = σ(x, t)``

```@example ONOFF
using ModelingToolkit, ModelingToolkitSampledData, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks
using JuliaSimCompiler
dt = 0.01
clock = Clock(dt)
z = ShiftIndex(clock)

qr(t) = 1sin(2t)  # reference position
qdr(t) = 2cos(2t) # reference velocity

function σ(x, t)
    q, qd = x
    e = q - qr(t)
    ė = qd - qdr(t)
    ė + e
end

@register_symbolic σ(x, t::Real)

@mtkmodel SuperTwistingSMC begin
    @structural_parameters begin
        z = ShiftIndex()
        Ts = SampleTime()
        nin
    end
    @components begin
        input = RealInput(; nin)
        output = RealOutput()
    end
    @parameters begin
        k = 1,    [description = "Control gain"]
        k2 = 1.1, [description = "Tuning parameter, often set to 1.1"]
    end
    @variables begin
        s(t) = 0.0, [description = "Sliding surface"]
        p(t) = 0.0
        y(t) = 0.0, [description = "Control signal output"]
        x(t) = 0.0
        xd(t) = 0.0
    end
    @equations begin
        s(z) ~ σ(input.u, t)
        p  ~ -√(k*abs(s))*sign(s)
        xd ~ -k2*k*sign(s)
        x(z) ~ x(z-1) + Ts*xd # Fwd Euler integration
        y ~ p + x
        output.u ~ y
    end
end


@mtkmodel ClosedLoopModel begin
    @components begin
        plant = SecondOrder(;)
        zoh = ZeroOrderHold()
        controller = SuperTwistingSMC(nin = 2, k=50)
    end
    @variables begin
        disturbance(t) = 0.0
    end 
    @equations begin
        controller.input.u ~ [Sample(dt)(plant.x), Sample(dt)(plant.xd)]
        disturbance ~ 2 + 2sin(3t) + sin(5t)
        connect(controller.output, zoh.input)
        zoh.output.u + disturbance ~ plant.input.u
    end
end
@named m = ClosedLoopModel()
m = complete(m)
ssys = structural_simplify(IRSystem(m))
prob = ODEProblem(ssys, [m.plant.x => -1, m.plant.xd => 0, m.controller.x(z-1) => 0], (0.0, 2π))
sol = solve(prob, Tsit5(), dtmax=0.01)
figy = plot(sol, idxs=[m.plant.x])
plot!(sol.t, qr.(sol.t), label="Reference")
figu = plot(sol, idxs=[m.zoh.y, m.disturbance], label=["Control signal" "Disturbance"])
plot(figy, figu, layout=(2,1))
```

The simulation indicates that the controller is able to track the reference signal despite the presence of the disturbance. The control signal exhibits a small degree of high-frequency chattering, a common characteristic of sliding-mode controllers.

