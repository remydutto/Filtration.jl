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
membrane_filtration_model(f₁::Function, f₂::Function, g::Function)
```
Note that respectively f_1, g and f₂ are tested to be L-functions and K-functions.

## Arguments
- f₁ : a L-function
- f₂ : a K-function
- g : a L-function

## Returns
- a membrane filtration model
"""
struct membrane_filtration_model
    # Parameters
    f₁::Function
    f₂::Function
    g::Function
    f₊::Function
    f₋::Function
    state_dynamic::Function
    cost_dynamic::Function
    
    function membrane_filtration_model(f₁::Function, f₂::Function, g::Function)
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
        f₊ = m -> 0.5 * (f₁(m) + f₂(m)) 
        f₋ = m -> 0.5 * (f₁(m) - f₂(m))
        state_dynamic = (m,u) -> f₋(m) + u * f₊(m)
        cost_dynamic = (m,u) -> u * g(m) 

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
function isLfunction(f::Function; start = 0, stop = 100, N = 100, ε = 10^-9)
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
function isKfunction(f::Function; start = 0, stop = 100, N = 100, ε = 10^-9)
    x = range(start, stop, N)
    y = f(x)
    !(ismonotonic(y,<)) && return false # Test of increasing value
    y < zero(y) && return false # Test of positive value on positive domain
    abs(f(0))> ε && return false # Test of lim = 0
    return true
end


function get_roots(model::membrane_filtration_model)
    @variables m
    Dₘ = Differential(m)

    g = model.g(m); f₊ = model.f₊(m); f₋ = model.f₋(m)

    Df₊ = Dₘ(f₊); Df₋ = Dₘ(f₋); Dg  = Dₘ(g)

    ν = g * (Df₋ * f₊ - f₋ * Df₊) + Dg * f₊ * f₋
    ν = simplify(expand_derivatives(ν))
    Dν = simplify(expand_derivatives(Dₘ(ν)))
    num, _ = Symbolics.arguments(Symbolics.value(ν))
    roots = symbolic_solve(num~0, m);
    roots = [Symbolics.symbolic_to_float(roots[i]) for i ∈ 1:length(roots)]
    ind = findall(x -> real(x)>0 && isreal(x), roots)
    roots = real(roots[ind])
    Dν = [substitute(Dν, m=>real(roots[i])) for i ∈ 1:length(roots)]
    return roots, Dν
end
