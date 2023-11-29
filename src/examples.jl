function solve_Poisson_equation(pos, ρ, Φ₀, Φₙ)
    # pos: positions
    # ρ:   space charge density in [C m^-3]
    # Φ₀:  potential at pos[1]
    # Φₙ:  potential at pos[1]

    N = length(pos)
    @assert length(pos) == length(ρ) # equal lengths
    
    dx = diff(pos)[1] 
    @assert all(δ ->  isapprox(δ,dx), diff(pos)) # equidistant grid

    # define Poisson equation in 1D with the boundary conditions Φ[1] = Φ₀, Φ[N] = Φₙ
    # central second derivative stencil
    a  = fill(1, N)
    b  = fill(-2, N)
    c  = fill(1, N)

    a[1] = 0   # has no 1st element
    c[N] = 0.0 # has no Nth element

    # define RHS
    ε₀ = 8.854188f-12
    d  = -dx^2/ε₀ .* ρ

    # left boundary 
    # equation: b[1] ⋅ x[1] + c[1] x[2] = d[1] 
    # Our boundary condition: x[1] = Φ₀
    # identify: c[1]=0, b[1]=1, d[1]=Φ₀
    b[1] = 1 
    c[1] = 0
    d[1] = Φ₀

    # right boundary 
    # equation: a[n] ⋅ x[n-1] + b[n] x[n] = d[n] 
    # Our boundary condition: x[n] = Φₙ
    # identify: a[n]=0, b[n]=1, d[n]=Φₙ
    b[N] = 1
    a[N] = 0
    d[N] = Φₙ

    # solve system
    Φ = Thomas_algorithm(a, b, c, d)
    return Φ  
end  