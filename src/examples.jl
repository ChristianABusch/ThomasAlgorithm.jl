function solve_Poisson_equation(x, f, g₀, gₙ)
    # solve the Poisson equation
    # Δu(x) = -f(x), with Diriclet boundary condition u(x[1])=g₀, u(x[end])=gₙ 

    N = length(x)
    @assert length(x) == length(f) # equal lengths
    
    dx = diff(x)[1] 
    @assert all(δ ->  isapprox(δ,dx), diff(x)) # equidistant grid

    # central second derivative stencil
    a  = fill(1,  N)
    b  = fill(-2, N)
    c  = fill(1,  N)

    a[1] = 0.0 # has no 1st element
    c[N] = 0.0 # has no Nth element

    # define RHS
    d  = -dx^2 .* f

    # left boundary 
    # equation: b[1] ⋅ u[1] + c[1] u[2] = d[1] 
    # Our boundary condition: u[1] = g₀
    # identify: c[1]=0, b[1]=1, d[1]=g₀
    b[1] = 1 
    c[1] = 0
    d[1] = g₀

    # right boundary 
    # equation: a[n] ⋅ u[n-1] + b[n] u[n] = d[n] 
    # Our boundary condition: u[n] = gₙ
    # identify: a[n]=0, b[n]=1, d[n]=gₙ
    b[N] = 1
    a[N] = 0
    d[N] = gₙ

    # solve system
    u = Thomas_algorithm(a, b, c, d)
    return u 
end  