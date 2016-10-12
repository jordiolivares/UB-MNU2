function jacobi(oldX, b)
    n = 10^6
    solution = zeros(n)
    solution[1] = (b[1] - oldX[3] - oldX[n-1])
    solution[2] = (b[2] - oldX[4] - oldX[n])
    for i = 3:n-2
        solution[i] = (b[i] - oldX[i-2] - oldX[i+2])
    end
    solution[n-1] = (b[n-1] - oldX[n-3] - oldX[1])
    solution[n] = (b[n] - oldX[n-2] - oldX[2])
    solution = solution * (1/3)
	solution
end

function gauss_seidel(oldX, b)
	sor(oldX, b, 1)
end

function sor_classe(oldX, b, ω)
	n = 10^6
	solution = zeros(n)
	solution[1] = oldX[1] + (ω/3) * (b[1] - oldX[3] - oldX[n-1])
	solution[2] = oldX[2] + (ω/3) * (b[2] - oldX[4] - oldX[n])
	for i = 3:n-2
		solution[i] = oldX[i] + (ω/3) * (b[i] - solution[i-2] - oldX[i+2])
	end
	solution[n-1] = oldX[n-1] + (ω/3) * (b[n-1] - solution[n-3] - solution[1])
	solution[n] = oldX[n] + (ω/3) * (b[n] - solution[n-2] - solution[2])
	return solution
end

function sor(oldX, b, ω)
	n = 10^6
	solution = zeros(n)
	solution[1] = (1-ω)oldX[1] + (ω/3) * (b[1] - oldX[3] - oldX[n-1])
	solution[2] = (1-ω)oldX[2] + (ω/3) * (b[2] - oldX[4] - oldX[n])
	for i = 3:n-2
		solution[i] = (1-ω)oldX[i] + (ω/3) * (b[i] - solution[i-2] - oldX[i+2])
	end
	solution[n-1] = (1-ω)oldX[n-1] + (ω/3) * (b[n-1] - solution[n-3] - solution[1])
	solution[n] = (1-ω)oldX[n] + (ω/3) * (b[n] - solution[n-2] - solution[2])
	return solution
end

function solve(convergence_cutoff, f::Function)
    n = 10^6
    oldX = randn(n)
    dif = Inf
    b = (1/n):(1/n):1
    β = 2/3 # variable temporal que fa que la norma es multipliqui per 1
    solution = zeros(n)
	iterations = 0
    while dif > convergence_cutoff
		# Canviar el mètode en aquí per sor_classe per a veure
		# com divergeix
        solution = f(oldX, b)
        dif = (β/(1-β)) * norm((solution - oldX), Inf)
        oldX = solution
		iterations += 1
    end
	println("Iterations: $iterations")
    return solution
end

sol_jacobi = solve(1/(10^12), jacobi)
sol_gauss = solve(1/(10^12), gauss_seidel)
sol_sor = solve(1/(10^12), (x, y) -> sor(x, y, 1.1))
println(sol_jacobi[1:5])
println(sol_jacobi[end-5:end])
println(sol_gauss[1:5])
println(sol_gauss[end-5:end])
println(sol_sor[1:5])
println(sol_sor[end-5:end])
