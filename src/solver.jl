function Thomas_algorithm!(a, b, c, d, x, c′, d′)
    N = length(b)

    # Forward sweep
    c′[1] = c[1]/b[1]
    for i in 2:N-1
        c′[i] = c[i]/(b[i] - a[i]* c′[i-1])
    end    

    d′[1] = d[1]/b[1]
    for i in 2:N
        d′[i] = (d[i] - a[i]* d′[i-1])/(b[i] - a[i]* c′[i-1])
    end  
    
    # Backward Sweep
    x[N] = d′[N]
    for i in N-1:-1:1
        x[i] = d′[i] - c′[i]*x[i+1]
    end
end

function Thomas_algorithm(a, b, c, d)
    N = length(b)

    x, c′, d′ = zeros(N), zeros(N), zeros(N)
    Thomas_algorithm!(a, b, c, d, x, c′, d′)
    return x
end