# set your working directory to the code
wd = "pn_kc/"
include(joinpath(wd,"modeling/learn.jl"))
# include("/usr/people/zhihaoz/scripts/miscellaneous/201006_pnkc_functions.jl")

function pick_gloms(act, num_act, nonact, num_nonact)
    return vcat(sample(act, num_act, replace=false),
    sample(nonact, num_nonact, replace=false))
end

include(joinpath(wd,"modeling/210625_pnkc_functions.jl"))

save_path = joinpath(wd,"modeling/")

# Figure7B-F
J = syn_conn_2k

Nkc = size(J,1)
Npn = NPN
date = "210627"


Ntrials = 10000
Nnetworks = 1000
epl = 0.4
P = 10
hp = Int(P/2)

conn_set = ["random_glom", "random_bouton", "random_claw","random_local","observed"]
group = "all_comm"
activated_dict = Dict("all_comm"=> all_comm, "core_comm"=>core_comm, "food"=>food)
t1 = activated_dict[group]
activated_gloms = findall(x -> x in t1, glom_seq)
nonactivated_set = setdiff(1:Npn, activated_gloms)
nlen = length(activated_gloms)

comm_fraction = [0, 4, 8, 12, 16, 19]

jn = length(conn_set)*length(comm_fraction)

r_mean = zeros(jn)
r_std = zeros(jn)
r_conn = String[]
r_commfract = Float32[]

ji = 1
for conn = conn_set
    for comm_pn = comm_fraction
        means = zeros(Nnetworks)
        for ni = 1:Nnetworks
            picked = pick_gloms(activated_gloms, comm_pn, nonactivated_set, nlen - comm_pn)
            if conn == "random_bouton"
                J = random_model(syn_conn_2k, sample_bouton)
            elseif conn == "random_glom"
                J = random_model(syn_conn_2k, sample_glom)
            elseif conn == "random_claw"
                J = random_claw(bi_conn_2k)
            elseif conn == "observed"
                J = vcat(syn_conn, syn_conn[rand(1:nkc_rd, 2000 - nkc_rd),:])
            else
                J = random_local(syn_conn, geom_tbl, nrows, row_idx)
            end
            means[ni] = mean(classify_s1(Npn, Nkc, nlen, picked, J, epl, P, hp, Ntrials))
            print(ni)
        end

        r_mean[ji] = mean(means)
        r_std[ji] = std(means)

        print("_mean")
        print(mean(means))
        print("_std")
        print(std(means))

        push!(r_conn, conn)
        push!(r_commfract, comm_pn)
        ji = ji + 1
        fname = join([date, conn, "5networks_10kTrials",epl, P, comm_pn],"-")
        npzwrite(string(save_path,"npz/",fname,".npz"), means)

    end
end

result_dict = Dict(
"means"=>r_mean,
"stds"=>r_std,
"conn"=> r_conn,
"NumComm"=>r_commfract
)
t11 = DataFrame(result_dict)

# CSV.write(string(save_path, fname, ".csv"), DataFrame(result_dict))


# Figure7G
#------------------------------------------------------

Nnetworks = 1000
Npn = NPN
Nkc = size(syn_conn_2k, 1)
Ntrials = 10000
P = 10
hp = Int(P/2)
date = "210624"
J = syn_conn_2k
act = "all_pn"
results = DataFrame([Float32, Float32, String, String, Float32, Int32],[:means,:stds,:conn,:act,:epl,:p], 300)
results[1:end,4] .= act
results[1:end,6] .= P
# epl_group = Array(0:0.1:1)
epl = 0.4
n = 1
act_group = ["all_pn"]
conn_group = ["random_glom", "random_bouton", "random_claw", "random_local","observed"]
for conn = conn_group
    means = zeros(Nnetworks)
    for ni = 1:Nnetworks
        if conn == "random_bouton"
            J = random_model(syn_conn_2k, sample_bouton)
        elseif conn == "random_glom"
            J = random_model(syn_conn_2k, sample_glom)
        elseif conn == "random_claw"
            J = random_claw(bi_conn_2k)
        elseif conn == "observed"
            J = vcat(syn_conn, syn_conn[rand(1:nkc_rd, 2000 - nkc_rd),:])
        else
            J = random_local(syn_conn, geom_tbl, nrows, row_idx)
        end
        means[ni] = mean(classify_s1(Npn, Nkc, Npn, 1:Npn, J, epl, P, hp, Ntrials))
    end
    fname = join([date, "conn", conn,
    "act", act,
    "epl", epl,
    "P", P],"-")
    npzwrite(string(save_path,"npz/",fname,".npz"), means)
    results[n,1] = mean(means)
    results[n,2] = std(means)
    print("_")
    print(std(means))
    print("_")
    print(mean(means))
    print("_")
    results[n,3] = conn
    results[n,5] = epl
    n = n + 1
    print(n)
end
fname = join([date, "act", act, "P", P],"-")
# CSV.write(string(save_path, "210626-allPNs_p10_noiseP4_5networks_1kTrials", ".csv"), results[1:5,1:end])



##------------------------------------------
# Supplemental Figure7H
J = syn_conn_2k

Nkc = size(J,1)
Npn = NPN
date = "210627"
epl = 0.4
P = 10
hp = Int(P/2)

conn_set = ["random_glom", "random_bouton", "random_claw","random_local","observed"]
group = "food"
activated_dict = Dict("all_comm"=> all_comm, "core_comm"=>core_comm, "food"=>food)
t1 = activated_dict[group]
activated_gloms = findall(x -> x in t1, glom_seq)
nonactivated_set = setdiff(1:Npn, activated_gloms)
nlen = length(activated_gloms)



# nlen=22
comm_fraction = [0, 22]

jn = length(conn_set)*length(comm_fraction)

r_mean = zeros(jn)
r_std = zeros(jn)
r_conn = String[]
r_commfract = Float32[]

ji = 1
for conn = conn_set
    for comm_pn = comm_fraction
        means = zeros(Nnetworks)
        for ni = 1:Nnetworks
            picked = pick_gloms(activated_gloms, comm_pn, nonactivated_set, nlen - comm_pn)
            if conn == "random_bouton"
                J = random_model(syn_conn_2k, sample_bouton)
            elseif conn == "random_glom"
                J = random_model(syn_conn_2k, sample_glom)
            elseif conn == "random_claw"
                J = random_claw(bi_conn_2k)
            elseif conn == "observed"
                J = vcat(syn_conn, syn_conn[rand(1:nkc_rd, 2000 - nkc_rd),:])
            else
                J = random_local(syn_conn, geom_tbl, nrows, row_idx)
            end
            means[ni] = mean(classify_s1(Npn, Nkc, nlen, picked, J, epl, P, hp, Ntrials))
            print(ni)
        end

        r_mean[ji] = mean(means)
        r_std[ji] = std(means)

        print("_mean")
        print(round(mean(means); digits=4))
        print("_std")
        print(round(std(means),digits=4))
        print("_")


        push!(r_conn, conn)
        push!(r_commfract, comm_pn)
        ji = ji + 1
        fname = join([date, conn, "10networks_1kTrials",epl, P, comm_pn],"-")
        npzwrite(string(save_path,"npz/",fname,".npz"), means)

    end
end

result_dict = Dict(
"means"=>r_mean,
"stds"=>r_std,
"conn"=> r_conn,
"NumComm"=>r_commfract
)
t11 = DataFrame(result_dict)

# CSV.write(string(save_path, fname, ".csv"), DataFrame(result_dict))


# SFigure7I
J = syn_conn_2k

Nkc = size(J,1)
Npn = NPN
date = "210627"



epl = 0.4
P = 10
hp = Int(P/2)

conn_set = ["random_glom", "random_bouton", "random_claw","random_local","observed"]
group = "core_comm"
activated_dict = Dict("all_comm"=> all_comm, "core_comm"=>core_comm, "food"=>food)
t1 = activated_dict[group]
activated_gloms = findall(x -> x in t1, glom_seq)
nonactivated_set = setdiff(1:Npn, activated_gloms)
nlen = length(activated_gloms)


# nlen=10
comm_fraction = [0, 10]

jn = length(conn_set)*length(comm_fraction)

r_mean = zeros(jn)
r_std = zeros(jn)
r_conn = String[]
r_commfract = Float32[]

ji = 1
for conn = conn_set
    for comm_pn = comm_fraction
        means = zeros(Nnetworks)
        for ni = 1:Nnetworks
            picked = pick_gloms(activated_gloms, comm_pn, nonactivated_set, nlen - comm_pn)
            if conn == "random_bouton"
                J = random_model(syn_conn_2k, sample_bouton)
            elseif conn == "random_glom"
                J = random_model(syn_conn_2k, sample_glom)
            elseif conn == "random_claw"
                J = random_claw(bi_conn_2k)
            elseif conn == "observed"
                J = vcat(syn_conn, syn_conn[rand(1:nkc_rd, 2000 - nkc_rd),:])
            else
                J = random_local(syn_conn, geom_tbl, nrows, row_idx)
            end
            means[ni] = mean(classify_s1(Npn, Nkc, nlen, picked, J, epl, P, hp, Ntrials))
            print(ni)
        end

        r_mean[ji] = mean(means)
        r_std[ji] = std(means)

        print("_mean")
        print(round(mean(means); digits=4))
        print("_std")
        print(round(std(means),digits=4))
        print("_")


        push!(r_conn, conn)
        push!(r_commfract, comm_pn)
        ji = ji + 1
        fname = join([date, conn, "5networks_1kTrials",epl, P, comm_pn],"-")
        npzwrite(string(save_path,"npz/",fname,".npz"), means)

    end
end

result_dict = Dict(
"means"=>r_mean,
"stds"=>r_std,
"conn"=> r_conn,
"NumComm"=>r_commfract
)
t11 = DataFrame(result_dict)

# CSV.write(string(save_path, fname, ".csv"), DataFrame(result_dict))


# SFig7J


J = syn_conn_2k

Nkc = size(J,1)
Npn = NPN
date = "210630"



epl = 0.4
P = 10
hp = Int(P/2)

conn_set = ["random_glom", "random_bouton", "random_claw","random_local","observed"]
group = "banana"
activated_dict = Dict("all_comm"=> all_comm, "core_comm"=>core_comm, "food"=>food, "banana"=>bas)
t1 = activated_dict[group]
activated_gloms = findall(x -> x in t1, glom_seq)
nonactivated_set = setdiff(1:Npn, activated_gloms)
nlen = length(activated_gloms)


comm_fraction = [0, 4]

jn = length(conn_set)*length(comm_fraction)

r_mean = zeros(jn)
r_std = zeros(jn)
r_conn = String[]
r_commfract = Float32[]

ji = 1
for conn = conn_set
    for comm_pn = comm_fraction
        means = zeros(Nnetworks)
        for ni = 1:Nnetworks
            picked = pick_gloms(activated_gloms, comm_pn, nonactivated_set, nlen - comm_pn)
            if conn == "random_bouton"
                J = random_model(syn_conn_2k, sample_bouton)
            elseif conn == "random_glom"
                J = random_model(syn_conn_2k, sample_glom)
            elseif conn == "random_claw"
                J = random_claw(bi_conn_2k)
            elseif conn == "observed"
                J = vcat(syn_conn, syn_conn[rand(1:nkc_rd, 2000 - nkc_rd),:])
            else
                J = random_local(syn_conn, geom_tbl, nrows, row_idx)
            end
            means[ni] = mean(classify_s1(Npn, Nkc, nlen, picked, J, epl, P, hp, Ntrials))
            print(ni)
        end

        r_mean[ji] = mean(means)
        r_std[ji] = std(means)

        print("_mean")
        print(round(mean(means); digits=4))
        print("_std")
        print(round(std(means),digits=4))
        print("_")


        push!(r_conn, conn)
        push!(r_commfract, comm_pn)
        ji = ji + 1
        fname = join([date, conn, "5networks_10kTrials",epl, P, comm_pn],"-")
        npzwrite(string(save_path,"npz/",fname,".npz"), means)

    end
end

result_dict = Dict(
"means"=>r_mean,
"stds"=>r_std,
"conn"=> r_conn,
"NumComm"=>r_commfract
)
t11 = DataFrame(result_dict)

CSV.write(string(save_path, fname, ".csv"), DataFrame(result_dict))
