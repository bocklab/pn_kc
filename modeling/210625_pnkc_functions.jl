
# set your working directory to the code
wd = "pn_kc/"
load_path = joinpath(wd,"modeling/data/")

using JLD
using DelimitedFiles
using Statistics
using CSV
using StatsBase
using DataFrames
using Random
using NPZ
using StatsBase

bi_conn = convert(Matrix, CSV.read(string(load_path,"200914-glom_kc_binary_conn.csv"), DataFrame, header=1)[:,2:end])

syn_list = convert(Matrix, CSV.read(string(load_path, "200914-ana_rd_synapses_per_claw.csv"), DataFrame, header=1)[:,2:end])

syn_conn = convert(Matrix, CSV.read(string(load_path, "200914-glom_kc_conn_w_syn.csv"), DataFrame, header=1)[:,2:end])

btn_cts = convert(Matrix, CSV.read(string(load_path,"200914-glom_kc_conn_bouton_count.csv"), DataFrame, header=1)[:,2:end])[:,2]

glom_seq = CSV.read(string(load_path,"200914-glom_seq.csv"), DataFrame, header=1)[:,2]

core_comm = ["DM2", "DP1m", "VM2", "DL2v", "DM3", "DM4", "DM1", "VM3", "VA2", "VA4"]

all_comm = vcat(core_comm, ["DP1l", "VA6", "DM6", "VC3l", "VM5d", "DL2d", "DM5", "VC4", "DC1"])

food = ["DM2", "DP1m", "VM2", "DL2v", "DM3", "DM4", "DM1", "VM3", "VA2", "VA4", "DM6", "DL2d", "VC3l", "DP1l", "VA6", "VM5d", "DC3", "VL2a", "VM1", "VC3m", "VM7d", "VM5v"]

bas = ["DM2","DM6","DM5","VM2"]

nkc_rd = size(syn_conn, 1)
syn_conn_2k = vcat(syn_conn, syn_conn[rand(1:nkc_rd, 2000 - nkc_rd),:])

nkc_rd = size(bi_conn, 1)
bi_conn_2k = vcat(bi_conn, bi_conn[rand(1:nkc_rd, 2000 - nkc_rd),:])

all_claws = sum(bi_conn, dims=2)
claws_2k = dropdims([Int(i) for i = sum(bi_conn_2k, dims=2)], dims=2)

num_pn = size(syn_conn_2k, 2)
num_kc = size(syn_conn_2k, 1)

t2 = StatsBase.proportionmap(dropdims(syn_list, dims=2))
const KS = collect(keys(t2))
const VALS = AnalyticWeights(collect(values(t2)))
const WTS = AnalyticWeights(btn_cts / sum(btn_cts))
const NPN = size(syn_conn_2k, 2)

function random_model(conn, sample_fun, remaining=claws_2k, ks=KS, vals=VALS, num_pn=NPN)
    # remember to drop claws singleton dimension
    # conn = conn_2k
    # claws = claws_2k
    # claws_2k, syn_list, num_pn
    rd_conn = zeros(size(conn))
    while sum(remaining) > 0
        picked_kcs = [ri[1] for ri in findall(!iszero, remaining)]
        picked_pns = rand(1:num_pn, length(picked_kcs))
        assign = CartesianIndex.(picked_kcs, sample_fun(length(picked_kcs)))
        rd_conn[assign] .+= sample(ks, vals, length(assign))
        remaining = @. (remaining > 0) * (remaining - 1)
    #    print(length(assign))
    end
    return rd_conn
end

sample_bouton(n) = sample(1:NPN, WTS, n)
sample_glom(n) = rand(1:NPN, n)

#random models
# random_model(conn, sample_bouton)
# random_model(conn, sample_glom)


# random claw null models
function random_claw(conn, ks=KS, vals=VALS)
    # conn = bi_conn_2k
    t9 = [findall(x->x>=i, conn) for i in 1:maximum(conn)];
    nnz = reduce(vcat, t9)
    new_nnz = CartesianIndex.(shuffle!([i[1] for i in nnz]), [i[2] for i in nnz])
    new_conn = zeros(size(conn))
    for i = new_nnz
    	new_conn[i] += sample(ks, vals)
    end
    return new_conn
end


function classify_s1(Npn, Nkc, nlen, picked, J, epl, P, hp, Ntrials)
    # Ntrialsper was set to 3
    # nlen: number of activated gloms
    # P, hp: number of targets and half number of targets
    # J: connectivity matrix
    # epl: noise levels
    # picked: index of picked activated gloms
    # epl_s = 1-epl
    epl_s = 1
    err = zeros(Ntrials)
    for ti = 1:Ntrials
        targets = ones(Int, P)
        targets[sample(1:P, hp, replace=false)] .= -1

        pats = zeros(Npn,P)
        pats[picked,:]= epl_s*randn(nlen,P)
        pats1 = pats + epl*randn(Npn,P)
        pats = (pats .> 0).*pats
        pats1 = (pats1 .> 0).*pats1
        mpats1 = J*pats1

        for mi = 1:Nkc
            mpats1[mi,:] .-= mean(mpats1[mi,:])
        end

        wm,bias= quadprogmargin(targets,mpats1)

        for jj = 1:3
            pats2 = pats + epl*randn(Npn,P)
            pats2 = (pats2 .> 0).*pats2
            mpats2 = J*pats2

            for mi = 1:Nkc
                mpats2[mi,:] .-= mean(mpats2[mi,:])
            end
            err[ti] += mean([sign(((mpats2)' * wm .+ bias)[i]) != targets[i] for i in 1:length(targets)])/3
        end
        err[err .> 0.99] .= 0.5
    end
    return err
end

t3 = CSV.read(string(load_path,"200926-RD_local_random_tbl.csv"), DataFrame, header=1)[:,2:end]

for rowi in 1:5
    t3[!, string("idx_",rowi)] = indexin(t3[!,rowi],glom_seq)
end
ids_for_claws = t3[!,:kc_ids]
kc_ids = unique(ids_for_claws)
row_idx = indexin(ids_for_claws, kc_ids)

geom_tbl = convert(Matrix, t3[!,["idx_1", "idx_2", "idx_3", "idx_4", "idx_5"]])
nrows = length(ids_for_claws)


function random_local(conn, geom_tbl, nrows, row_idx, ks=KS, vals=VALS)
    # conn: the small subset of KCs manually traced only
    col_idx = geom_tbl[CartesianIndex.(1:nrows, sample(1:5,nrows))]
    new_conn = zeros(size(conn))
    for ci in CartesianIndex.(row_idx, col_idx)
        new_conn[ci] += sample(ks,vals)
    end
    nkc_rd = size(conn, 1)
    new_conn_2k = vcat(new_conn, new_conn[rand(1:nkc_rd, 2000 - nkc_rd),:])
    return new_conn_2k
end
