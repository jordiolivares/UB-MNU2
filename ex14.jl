function jacobi(x̄, b)
	n = 10^6
	x = zeros(n)
	x[1] = (b[1] - x̄[3] - x̄[n-1])
	x[2] = (b[2] - x̄[4] - x̄[n])
	for i = 3:n-2
		x[i] = (b[i] - x̄[i-2] - x̄[i+2])
	end
	x[n-1] = (b[n-1] - x̄[n-3] - x̄[1])
	x[n] = (b[n] - x̄[n-2] - x̄[2])
	x = x * (1/3)
	x
end

function gauss_seidel(x̄, b)
	sor(x̄, b, 1)
end

function sor_classe(x̄, b, ω)
	n = 10^6
	x = zeros(n)
	x[1] = x̄[1] + (ω/3) * (b[1] - x̄[3] - x̄[n-1])
	x[2] = x̄[2] + (ω/3) * (b[2] - x̄[4] - x̄[n])
	for i = 3:n-2
		x[i] = x̄[i] + (ω/3) * (b[i] - x[i-2] - x̄[i+2])
	end
	x[n-1] = x̄[n-1] + (ω/3) * (b[n-1] - x[n-3] - x[1])
	x[n] = x̄[n] + (ω/3) * (b[n] - x[n-2] - x[2])
	x
end

function sor(x̄, b, ω)
	n = 10^6
	x = zeros(n)
	x[1] = (1-ω)x̄[1] + (ω/3) * (b[1] - x̄[3] - x̄[n-1])
	x[2] = (1-ω)x̄[2] + (ω/3) * (b[2] - x̄[4] - x̄[n])
	for i = 3:n-2
		x[i] = (1-ω)x̄[i] + (ω/3) * (b[i] - x[i-2] - x̄[i+2])
	end
	x[n-1] = (1-ω)x̄[n-1] + (ω/3) * (b[n-1] - x[n-3] - x[1])
	x[n] = (1-ω)x̄[n] + (ω/3) * (b[n] - x[n-2] - x[2])
	x
end

function solve(convergence_cutoff, f::Function)
	n = 10^6
	x̄ = randn(n)
	dif = Inf
	b = (1/n):(1/n):1
	β = 2/3 # norma de la matriu
	x = zeros(n)
	iterations = 0
	while dif > convergence_cutoff
		x = f(x̄, b)
		dif = (β/(1-β)) * norm((x - x̄), Inf)
		x̄ = x
		iterations += 1
	end
	println("Iterations: $iterations")
	return x
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
println("SOR (ω = 0.8):")
println()
map(x -> @printf("%.12f\n", x), sol_sor[1:10])
println()
