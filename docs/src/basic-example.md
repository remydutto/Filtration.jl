# Membrane Filtration : Direct solve

Let us consider the membrane filtration process where the goal is to maximize the water net production 

```math
    \max_{x(\cdot), u(\cdot)} \int_{0}^{T} u(t) g\big( m(t) \big) \mathrm dt, 
```

where the control $u(\cdot) \in \mathrm L^\infty([0,T], \mathbb R)$ corresponds to the filtration mode (1 during filtration and -1 during backwash) and the state $m (\cdot) \in \mathrm{AC}([0, T], \mathbb R)$ is the mass of the cake layer formed during the water filtration. The dynamic of this state is given by 

```math
    \dot m(t) = \frac{1 + u(t)}{2} f_1(m(t)) - \frac{1 - u(t)}{2} f_2(m(t))
```
with the initial condition $m(0) = m_0$. 

Based on [[Kalboussi et al., 2018](https://doi.org/10.1016/j.ifacol.2017.08.1554)], we assume that $f_1$ and $g$ are smooth $\mathcal L$-function, and that $f_2$ is a smooth $\mathcal K$ function. 

First, we need to import the [Filtration.jl](https://remydutto.github.io/Filtration.jl) to define this Optimal Control Problem, to give it to [OptimalControl.jl](https://control-toolbox.org/OptimalControl.jl) package and [NLPModelsIpopt.jl](jso.dev/NLPModelsIpopt.jl) to solve it. 
We also need to import the [Plots.jl](https://docs.juliaplots.org) package to plot the solution.

```@example main
using Filtration
using CTBase
using OptimalControl
using Plots
using MadNLP
```

## Model definition 

In this example, these functions are defined by using the [[Benyahia et al.](https://www.researchgate.net/publication/272506325_A_simple_model_of_anaerobic_membrane_bioreactor_for_control_design_coupling_the_AM2b_model_with_a_simple_membrane_fouling_dynamics)] model : 

```math
f_1(m) = \frac{b}{e+m}, \quad f_2(m) = am, \quad f_3(m) = \frac{1}{e+m}, 
```
where $a$, $b$ and $e$ are positive numbers. The problem is simply defined, thanks to `MembraneFiltrationModel`.

```@example main
# Benyahia & al. model
a = 1; b = 1; e = 1;
f₁ = m ->  b ./ (e .+ m)
f₂ = m -> a .* m
g  = m -> 1 ./ (e .+ m)
model = MembraneFiltrationModel(f₁, f₂, g);
```

## Verification of hypotheses

With this model, we know that the function $\psi$ defined by 

```math
\psi(m) = g(m)\big[f_{-}'(m) f_{+}(m) - f_{-}(m) f_{+}'(m)\big] + g'(m) f_+(m) f_-(m)
```

have an unique positive root, denoted $\bar m$, where function $f_+$ and $f_-$ are given by 

```math
f_+(m) = \frac{f_1(m) + f_2(m)}{2} \quad f_-(m) = \frac{f_1(m) -f_2(m)}{2} \cdot
```

For more information about why this property is needed, please refer to [[Kalboussi et al., 2018](https://doi.org/10.1016/j.ifacol.2017.08.1554)]. However, we can take a look on this function $\psi$. 


```@example main
ψ = get_psi(model)
root = get_root(model)
plot(range(0,10, 100), ψ, label = "ψ")
scatter!([root], [ψ(root)], label = "zero of ψ")
```

## Direct Solve

We can solve the optimal control problem by derect method, thanks to the [OptimalControl.jl](https://control-toolbox.org/OptimalControl.jl) package.

```@example main
t0 = 0; m0 = 1; tf = 10;                # initial and final time and state
@def ocp begin                          # problem definition
    t ∈ [ t0, tf ], time
    m ∈ R, state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    m(t0) == m0
    ṁ(t) == model.state_dynamic(m(t),u(t))
    ∫( model.cost_dynamic(m(t),u(t))) → max
end
sol = OptimalControl.solve(ocp, :madnlp)
plot(sol)
```