include("learn.jl")
using CSV
using StatsBase
using DataFrames
using Statistics

# the observed connectivity, with first 16 PNs being the community PNs
# (from 190421-process_pnkc_connectivity.py)
ro_conn = convert(Matrix, CSV.read("data/190319-440RDKC_PNfirst16comm_last22blue.csv", header=1)[:,2:end])

Ntrials = 100
Ntrialsper = 3
Nkc = size(ro_conn,1)
Npn = size(ro_conn,2)
f = 0.05
P = 10
eps = 0.2

# as of now, running the following takes quite a while
# number of PNs as part of the 16 with acivity patterns
num_comm_pns = [0, 2, 4, 6, 8, 10, 12, 14, 16]

result_dict = Dict()
J = ro_conn

num_exp = 100
for comm_pn = num_comm_pns
	mean_set = Any[]
	std_set = Any[]
	sqrt_set = Any[]
	for t2 = 1:num_exp
        picked = vcat(sample(1:16, comm_pn, replace=false),
        sample(17:Npn, 16 - comm_pn, replace=false))
        err = zeros(Ntrials)
        for ti = 1:Ntrials
            targets = rand(-1:2:1,P)

            pats = zeros(Npn,P)
            pats[picked,:]= randn(16,P)

            pats1 = pats + eps*randn(Npn,P)

            pats = (pats .> 0).*pats
            pats1 = (pats1 .> 0).*pats1

            h1 = J*pats1

            thresh = zeros(Nkc)
            mpats1 = zeros(Nkc,P)
            for ci = 1:Nkc
                allh = sort(h1[ci,:])
                thresh[ci] = (allh[floor(Int,(1-f)*length(allh))] + allh[floor(Int,(1-f)*length(allh))+1])/2
                mpats1[ci,:] = (h1[ci,:] .> thresh[ci]) .* (h1[ci,:] .- thresh[ci])
            end

            for mi = 1:Nkc
                mpats1[mi,:] .-= mean(mpats1[mi,:])
            end

            wm,bias= quadprogmargin(targets,mpats1)

            for jj = 1:Ntrialsper
                pats2 = pats + eps*randn(Npn,P)
                pats2 = (pats2 .> 0).*pats2
                h2 = J*pats2
                mpats2 = zeros(Nkc,P)

                for ci = 1:Nkc
                    mpats2[ci,:] = (h2[ci,:] .> thresh[ci]) .* (h2[ci,:] .- thresh[ci])
                end

                for mi = 1:Nkc
                    mpats2[mi,:] .-= mean(mpats2[mi,:])
                end

                err[ti] += mean([sign(((mpats2)' * wm .+ bias)[i]) != targets[i] for i in 1:length(targets)])/Ntrialsper
            end
        end
        err[err .> 0.99] .= 0.5
        push!(mean_set, mean(err))
        push!(std_set, std(err))
        push!(sqrt_set, std(err)/sqrt(length(err)))
        print(t2)
	end
	result_dict[string(comm_pn)] = copy(mean_set)
end

using DataFrames
df = DataFrame(result_dict)

# CSV.write("190422-comm_fract_vs_performance.csv", df)
