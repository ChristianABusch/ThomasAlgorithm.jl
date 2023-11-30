# ThomasAlgorithm

This package implements the tridiagonal matrix algorithm, also known as the Thomas algorithm, which solves a system of linear equations with a tridiagonal matrix. Variables are defined as described here: [Wikipedia](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)

## API

```julia
using ThomasAlgorithm
a = [0, 3, 2] # lower diagonal
b = [1, 1, 1] # diagonal
c = [1, 1, 0] # upper diagonal
d = [2, 0, 4] # right hand side

# solve for x:
x = Thomas_algorithm(a,b,c,d)
```

## Solving the Poisson equation
A function for solving the Poisson equation 
$$\frac{d^2 u}{d x^2} = - f(x)$$,
using the Thomas algorithm, is also implemented. 

### Example
The solution to the problem
$$\frac{d^2 u}{d x^2} = - f(x),  \quad x \in [-10, 10], \quad f(x)=\cos(\frac{\pi}{2L}x)$$
with the Diriclet boundary condition
$$u(-10)=-10, u(10) = 20$$
can be found as follows:

```julia
    N = 50
    L = 10.0
    x = collect(LinRange(-L,L,N))
    u₀, uₙ = -10,20
    f   = @. cos(π/(2*L)*x)
	
    u= ThomasAlgorithm.solve_Poisson_equation(x, f, u₀, uₙ)
```
