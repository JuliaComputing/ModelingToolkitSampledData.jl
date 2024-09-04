# On-off controller
A common form of controller is one which can output two distinct values only, such as on/off, or 1/-1. Below, we demonstrate the usage of such a controller, the [`DiscreteOnOffController`](@ref), and let it control an unstable first-order system. The system is given by the equation
```math
τ \dot x = x + 10u
```
where ``τ = 100`` is a time constant. The controller ``u = f(x)`` is an on/off controller with a hysteresis band of 0.3 and possible output values -1 and 1. The controller runs with an update frequency of 1 Hz.


```@example ONOFF
cl = Clock(1)
z = ShiftIndex(cl)
@mtkmodel OnOffModel begin
    @components begin
        onoff = DiscreteOnOffController(; b = 0.3, z, bool=false)
        c = Constant(k=1)
    end
    @variables begin
        x(t) = 0 
    end
    @equations begin
        onoff.u ~ Sample(cl)(x)
        connect(c.output, onoff.reference)
        100D(x) ~ x + 10Hold(onoff.y)
    end
end
@named m = OnOffModel()
m = complete(m)
ssys = structural_simplify(IRSystem(m))
prob = ODEProblem(ssys, [m.onoff.y(z-1) => 0], (0.0, 40.0))
sol = solve(prob, Tsit5(), dtmax=0.1)
plot(sol, idxs=[m.x, m.onoff.y], title="On-off control of an unstable first-order system")
```

The parameter `b` controls the width of the hysteresis band, and the structural parameter `bool = false` indicates that this is a 1/-1 controller and not a 0/1 controller. The controller additionally takes a parameter `k` which is the gain of the controller, here we used the default `k = 1`.