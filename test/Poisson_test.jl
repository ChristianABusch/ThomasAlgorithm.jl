function Poisson_equation_test(atol)
   # define physical grid
    N = 50
    L = 10.0
    x = collect(LinRange(-L,L,N))

    # Problem: 
    # -d²u/dx² = f(x),   x ∈ [-L,L]
    # f = f₀ cos(π/(2L)⋅x) 
    # u(-L)=-10, u(L)=20
    u₀, uₙ = -10,20

    # Numerical solution
    f   = @. cos(π/(2*L)*x)
    sol_numerical = ThomasAlgorithm.solve_Poisson_equation(x, f, u₀, uₙ)

    # analytic solution
    c₂ = -(u₀ + uₙ)/2
    c₁ = -(uₙ + c₂)/L
    sol_analytical = @. (2*L/π)^2 * cos(π/(2*L)*x) - c₁ * x - c₂

    return all(abs.(sol_analytical .-  sol_numerical) .< atol)
end
