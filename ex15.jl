# Calculates x^t * A
function transpA(x)
	n = 10^6
	vector = zeros(n)
	vector[1] = 3x[1] + x[3] + x[n-1]
	vector[2] = 3x[2] + x[4] + x[n]
	for i = 3:n-2
		vector[i] = 3x[i] + x[i+2] + x[i-2]
	end
	vector[n-1] = 3x[n-1] + x[n-3] + x[1]
	vector[n] = 3x[n] + x[n-2] + x[2]
	vector
end

function ∇Q(x, b)
	vector = transpA(x)
	vector - b
end

function steepestDescent(oldX, b, ω)
	p_k = ∇Q(oldX, b)
	a_k = dot(p_k, p_k)/dot(transpA(p_k), p_k)
	x = oldX - ω * a_k * p_k
	x
end

function steepestDescent(oldX, b)
	steepestDescent(oldX, b, 1)
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
	return x
end

sol = solve(1/10^12, (x,y) -> steepestDescent(x, y, 0.8))
println("Steepest Descent (ω = 0.8):")
println()
map(x -> @printf("%.12f\n", x), sol[1:10])
println()
