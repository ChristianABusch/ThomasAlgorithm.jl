function Poisson_equation_test(atol)
    # define physical grid
    N = 50
    L = 0.01
    x = collect(LinRange(-L,L,N))

    # Problem: 
    # d²Φ/dx² = -ρ/ε₀,   x ∈ [-L,L]
    # ρ = e ⋅ n₀ cos(π/(2L)⋅x) 
    # Φ(-L)=-10, Φ(L)=20
    n₀     = 1e14
    q,  ε₀ = 1.602e-19, 8.854e-12
    Φ₀, Φₙ = -10,20

    # Numerical solution
    ρ   = @. q * n₀ * cos(π/(2*L)*x)
    sol_numerical = ThomasAlgorithm.solve_Poisson_equation(x, ρ, Φ₀, Φₙ)

    # analytic solution
    c₂ = -(Φ₀ + Φₙ)/2
    c₁ = -(Φₙ + c₂)/L
    sol_analytical = @. q * n₀/ε₀  * (2*L/π)^2 * cos(π/(2*L)*x) - c₁ * x - c₂

    return all(abs.(sol_analytical .-  sol_numerical) .< atol)
end
