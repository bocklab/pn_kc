# originally from Ashok Litwin-Kumar (Eichler et al. 2017)
# with slight modification (refer to larva/learn.jl)

using Ipopt,MathProgBase,LinearAlgebra
#return maximal margin solution
function quadprogmargin(labels,pats)
	N = size(pats,1)
	P = length(labels)

#	Q = diagm(ones(N+1))
	Q = convert(Array, Diagonal(ones(N+1)))
	Q[N+1,N+1] = 0

	c = zeros(N+1)

	b = ones(P)

	A = zeros(N+1,P)

	for ni = 1:N
		for qi = 1:P
			A[ni,qi] = pats[ni,qi]*labels[qi]
		end
	end
	A[N+1,:] = labels

	sol = quadprog(c,Q,A','>',b,-Inf,Inf,IpoptSolver(print_level=0,max_iter=400))
	if length(sol.sol) > 0
		w = sol.sol[1:N]
		bias = sol.sol[N+1]
		return w/norm(w),bias/norm(w)
	else
		return zeros(N),0
	end
end
