# Membrane Filtration Model

"""
$(TYPEDSIGNATURES)

A membrane filtration model.

# Contains:
- a L-function f₁
- a K-function f₂
- a L-function g
- the function f₊
- the function f₋
- a function which describes the dynamic of the state (x,u) -> f(x,u)
- a function which describes the dynamic of the cost (x,u) -> f⁰(x,u)

# Constructor
A membrane filtration model can be constructed using the following code:
```julia
MembraneFiltrationModel(f₁::Function, f₂::Function, g::Function)
```
Note that respectively f_1, g and f₂ are tested to be L-functions and K-functions.

## Arguments
- f₁ : a L-function
- f₂ : a K-function
- g : a L-function

## Returns
- a MembraneFiltrationModel
"""
struct MembraneFiltrationModel
    # Parameters
    f₁::Function
    f₂::Function
    g::Function
    f₊::Function
    f₋::Function
    state_dynamic::Function
    cost_dynamic::Function
    
    function MembraneFiltrationModel(f₁::Function, f₂::Function, g::Function)
        if !(isLfunction(f₁))
            display(L"Please verify that $f_1$ is a smooth $\mathcal L$-function : decreasing with $\lim_{x \to \infty} f_1(x) = 0$.")
            display(plot(f₁, label = "f₁", size = (500, 300)))
            error("Wrong definition of inputs functions")
        end
        if !(isLfunction(g))
            display(L"Please verify that $g$ is a smooth $\mathcal L$-function : decreasing with $\lim_{x \to \infty} g(x) = 0$.")
            display(plot(g, label = "g", size = (500, 300)))
            error("Wrong definition of inputs functions")
        end
        if !(isKfunction(f₂))
            display(L"Please verify that $f₂$ is a smooth $\mathcal K$-function : increasing with $g(0) = 0$.")
            display(plot(f₂, label = "f₂", size = (500, 300)))
            error("Wrong definition of inputs functions")
        end
        f₊(m) = 0.5 * (f₁(m) + f₂(m)) 
        f₋(m) = 0.5 * (f₁(m) - f₂(m))
        state_dynamic(m,u) = f₋(m) + u * f₊(m)
        cost_dynamic(m,u) = u * g(m) 

        model = new(
            f₁,
            f₂,
            g,
            f₊,
            f₋,
            state_dynamic,
            cost_dynamic
        )
        return model
    end
end

"""
$(TYPEDSIGNATURES)

Check if a vector is monotonic, with respect to a given comparison operator.

# Arguments
- V : a vector
- cmp : a comparison operator, initialized to >

# Returns
- a boolean

"""
function ismonotonic(V::AbstractVector, cmp = >)
    current = V[begin]
    for i ∈ 1:length(V)-1
        newval = V[i+1]
        !cmp(current, newval) && return false
        current = newval
    end
    return true
end

"""
$(TYPEDSIGNATURES)

Check if a function is a L-function.

# Arguments
- f : a function
- start : the start of the domain, initialized to 0
- stop : the end of the domain, initialized to 100
- N : the number of points, initialized to 100
- ε : the precision, initialized to 10^-9

# Returns
- a boolean

"""
function isLfunction(f::Function; start::Int = 0, stop::Int = 100, N::Int = 100, ε::Real = 10^-9)
    x = range(start, stop, N)
    y = f(x)
    !(ismonotonic(y,>)) && return false # Test of decreasing value   
    y < zero(y) && return false # Test of positive value on positive domain
    f(10^12)> ε && return false # Test of lim = 0
    return true
end

"""
$(TYPEDSIGNATURES)

Check if a function is a K-function.

# Arguments
- f : a function
- start : the start of the domain, initialized to 0
- stop : the end of the domain, initialized to 100
- N : the number of points, initialized to 100
- ε : the precision, initialized to 10^-9

# Returns
- a boolean

"""
function isKfunction(f::Function; start::Int = 0, stop::Int = 100, N::Int = 100, ε::Real = 10^-9)
    x = range(start, stop, N)
    y = f(x)
    !(ismonotonic(y,<)) && return false # Test of increasing value
    y < zero(y) && return false # Test of positive value on positive domain
    abs(f(0))> ε && return false # Test of lim = 0
    return true
end

"""
$(TYPEDSIGNATURES)

Construct the function ψ by using ForwardDiff.

# Arguments
- model : a MembraneFiltrationModel

# Returns
- ψ : the function ψ
"""
function get_ψ(model::MembraneFiltrationModel)
    f₊, f₋, g = model.f₊, model.f₋, model.g
    df₊(m) = ForwardDiff.derivative(f₊, m)
    df₋(m) = ForwardDiff.derivative(f₋, m)
    dg(m) = ForwardDiff.derivative(g, m)
    ψ(m) = g(m) * (df₋(m) * f₊(m) - f₋(m) * df₊(m)) + dg(m) * f₊(m) * f₋(m)
    return ψ
end

"""
$(TYPEDSIGNATURES)

Compute the positive roots of the function ψ and its derivative by using Symbolics. 
This function works only if the function ψ is an algebraic fraction.

!!! note
    This function may be long to compute due to the use of Symbolics. 
    However, if it compute, it assure that all positive roots of ψ are given.
    If ψ have only one positive root, use get_roots instead.

# Arguments
- model : a MembraneFiltrationModel

# Returns
- roots : the positive roots of ψ
- Dψ : the derivative of ψ at the roots

"""
function get_roots_symbolic_algebraic_fraction(model::MembraneFiltrationModel)
    @variables m
    Dₘ = Differential(m)
    g = model.g(m); f₊ = model.f₊(m); f₋ = model.f₋(m)
    Df₊ = Dₘ(f₊); Df₋ = Dₘ(f₋); Dg  = Dₘ(g)
    ψ = g * (Df₋ * f₊ - f₋ * Df₊) + Dg * f₊ * f₋
    ψ = simplify(expand_derivatives(ψ))
    Dψ = simplify(expand_derivatives(Dₘ(ψ)))
    num, _ = Symbolics.arguments(Symbolics.value(ψ))
    roots = symbolic_solve(num~0, m);
    roots = [Symbolics.symbolic_to_float(roots[i]) for i ∈ 1:length(roots)]
    ind = findall(x -> real(x)>0 && isreal(x), roots)
    roots = real(roots[ind])
    Dψ = [substitute(Dψ, m=>real(roots[i])) for i ∈ 1:length(roots)]
    return roots, Dψ
end

"""
$(TYPEDSIGNATURES)

Compute a root of the function ψ by using ForwardDiff.

!!! note
    This function return a root of ψ. 
    If ψ have multiple roots, it may not return the desired one.

# Arguments
- model : a MembraneFiltrationModel
- x₀ : the initial guess, initialized to 0.5

# Returns
- roots : a root of ψ

"""
function get_root(model::MembraneFiltrationModel, x₀::Real = 0.5)
    ψ = get_ψ(model)
    root = find_zero(ψ, x₀)
    return root
end
