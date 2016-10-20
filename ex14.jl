function jacobi(oldX, b)
	n = 10^6
	x = zeros(n)
	x[1] = (b[1] - oldX[3] - oldX[n-1])
	x[2] = (b[2] - oldX[4] - oldX[n])
	for i = 3:n-2
		x[i] = (b[i] - oldX[i-2] - oldX[i+2])
	end
	x[n-1] = (b[n-1] - oldX[n-3] - oldX[1])
	x[n] = (b[n] - oldX[n-2] - oldX[2])
	x = x * (1/3)
	x
end

function gauss_seidel(oldX, b)
	sor(oldX, b, 1)
end

function sor_classe(oldX, b, ω)
	n = 10^6
	x = zeros(n)
	x[1] = oldX[1] + (ω/3) * (b[1] - oldX[3] - oldX[n-1])
	x[2] = oldX[2] + (ω/3) * (b[2] - oldX[4] - oldX[n])
	for i = 3:n-2
		x[i] = oldX[i] + (ω/3) * (b[i] - x[i-2] - oldX[i+2])
	end
	x[n-1] = oldX[n-1] + (ω/3) * (b[n-1] - x[n-3] - x[1])
	x[n] = oldX[n] + (ω/3) * (b[n] - x[n-2] - x[2])
	x
end

function sor(oldX, b, ω)
	n = 10^6
	x = zeros(n)
	x[1] = (1-ω)oldX[1] + (ω/3) * (b[1] - oldX[3] - oldX[n-1])
	x[2] = (1-ω)oldX[2] + (ω/3) * (b[2] - oldX[4] - oldX[n])
	for i = 3:n-2
		x[i] = (1-ω)oldX[i] + (ω/3) * (b[i] - x[i-2] - oldX[i+2])
	end
	x[n-1] = (1-ω)oldX[n-1] + (ω/3) * (b[n-1] - x[n-3] - x[1])
	x[n] = (1-ω)oldX[n] + (ω/3) * (b[n] - x[n-2] - x[2])
	x
end

function solve(convergence_cutoff, f::Function)
	n = 10^6
	solution = randn(n)
	dif = Inf
	b = (1/n):(1/n):1
	β = 2/3 # norma de la matriu
	x = zeros(n)
	iterations = 0
	while dif > convergence_cutoff
		x = f(solution, b)
		dif = (β/(1-β)) * norm((x - solution), Inf)
		solution = x
		iterations += 1
	end
	println("Iterations: $iterations")
	x
end

sol_jacobi = solve(1/(10^12), jacobi)
sol_gauss = solve(1/(10^12), gauss_seidel)
sol_sor = solve(1/(10^12), (x, y) -> sor(x, y, 1.1))
println("Jacobi:")
println()
map(x -> @printf("%.12f\n", x), sol_jacobi[1:10])
println()
println("Gauss-Seidel:")
println()
map(x -> @printf("%.12f\n", x), sol_gauss[1:10])
println()
println("SOR (ω = 1.1):")
println()
map(x -> @printf("%.12f\n", x), sol_sor[1:10])
println()
