
# This is modiried and extended from Eichler et al. code provided by Ashok Litwin-Kumar
# updated to use the newest Julia version (v1.1)

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

# coding level
f = 0.05

# number of odors
P = 10

# standard deviation of Gaussian noise
eps = 0.2

# randomized connectivity based on number of PNs
rd_conn = zeros(size(ro_conn))
for i = 1:size(rd_conn,1)
	t1 = ro_conn[i,:]
	nz = t1[t1.!=0]
	rd_conn[i, sample(1:Npn, length(nz), replace=false)] = nz
end

# observed and random connectivity matrix
J_set = Dict("observed" => ro_conn, "random" => rd_conn)

results = Any[]
for key = ["observed","random"]

	# for observed network, compare between community and non-community
	# for random network, only one (community or non-community is equivilent
	# since connectivity is randomized)
	if key == "observed"
		compare_comm = [true, false]
	else
		compare_comm = true
	end

	J = J_set[key]

	for compare = compare_comm
		err = zeros(Ntrials)
		for ti = 1:Ntrials
			targets = rand(-1:2:1,P)

			if compare

				# in the connectivity matrix, I put the communiy PNs in the first 16 columns, so when activity patterns are given from the first 16 PNs, they are created in the community PNs
				pats = vcat(randn(16,P), zeros(Npn-16, P))
			else

				# otherwise, select 16 PNs from the non-community ones randomly
				pats = zeros(Npn,P)
				pats[rand(17:Npn,16),:]= randn(16,P)

			end
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

			# use the labels (targets) and activity patterns (mpats1) to train the classifier
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

				# compare predictions by the classifier and same input patterns with different noise realizations
				err[ti] += mean([sign(((mpats2)' * wm .+ bias)[i]) != targets[i] for i in 1:length(targets)])/Ntrialsper

			end
		end

		# it is a binary categorization task, so chance is 0.5
		err[err .> 0.99] .= 0.5

		push!(results,err)

        print(mean(err))
	end
end
results

df = DataFrame(ob_comm=results[1], ob_noncomm=results[2], rd=results[3])

# conclusion (1a): mean(ob_comm) < mean(ob_non_comm)
# conclusion (1b): mean(ob_comm) < mean(rd), only slightly, not robust

# Notes:
# 1) The conclusions are:
#    a) the observed PN-KC connectivity leads to discrimination performance of
#    community input patterns better than that of non-community
#    input patterns;
#    b) compared to random connectivity, the observed PN-KC connectivity
#    leads to slightly better performance in discriminating community input
#    patterns.

# 2) Conclusion (a) is robust. Namely, the difference between mean(ob_comm) and
#    mean(ob_noncomm). However, (b) is NOT robust and only slightly better.

# 3) The comparison may not be equivilent since
#	for ob_comm, it is averaging over the same set of PNs and
#   for ob_noncomm, it is averaging over different set of randomly select PNs.
#   However, I do have a different approach to address this issue
#   (190422_varying_comm_fractions.jl).
