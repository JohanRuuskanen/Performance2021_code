
function extractParameters(suffix::String, subCN::Array{String, 2})
    matfile = matopen(joinpath(datapath, "params_$suffix.mat"))

    Queues = String.(read(matfile, "Queues")[:])
    Classes = String.(read(matfile, "Classes")[:])
    Disc = uppercase.(read(matfile, "Disc")[:])
    Q0 = Int64.(read(matfile, "Q0"))
    Servers = read(matfile, "K")[:]
    Pa = read(matfile, "Pa")
    ph_params = read(matfile, "E")

    @assert all((s -> any(occursin.(s, ["INF", "PS", "EXT"]))).(Disc))

    # Amount of non source/sink Queues
    Q = length(Queues) - sum(Disc .== "EXT")

    # Get transition matrix
    P_all = P_NCP_splat(classQueueP(Pa, Queues, Classes))

    idx_trans, idx_sources, idx_sinks = getUsedNodes(P_all, Q0)

    P = P_all[idx_trans, idx_trans]
    E = P_all[idx_sources, idx_trans]

    # Get the servers
    K = Servers[Disc .!= "EXT"]

    # Start values
    X0 = Float64.((Q0'[:])[idx_trans])

    # Classes per queue
    which_class = [idx_trans[(k-1)*length(Classes)+1:k*length(Classes)] for k = 1:length(Queues)]
    C = sum.(which_class)[Disc .!= "EXT"]
    C_idx = [sum(C[1:k-1])+1:sum(C[1:k]) for k in 1:Q]
    C_set = (wc -> Classes[wc]).(which_class)[Disc .!= "EXT"]
    @assert length(C) == Q

    # Extract PH dists
    ph, ph_arrivals = getPHFromMatrix(ph_params, sum(C), (Disc[1] == "EXT"))

    # Read phase distribution
    S, M, Mc, N = getPhasesPos(ph, Q, C)

    # Extract the CN set indexes and transition matrix
    subCNa = zeros(Int64, size(subCN, 1))
    for k = 1:length(subCNa)
        iq = findfirst(subCN[k, 2] .== Queues[Disc .!= "EXT"])
        ic = findfirst(subCN[k, 1] .== C_set[iq])
        subCNa[k] = sum(C[1:iq-1]) + ic
    end
    subCNaC = filter(x -> x ∉ subCNa, 1:sum(C))

    # Create parameter struct
    Params = ParameterStruct(Queues, Classes, Disc, Q, C, C_idx, C_set, K, 
        subCN, subCNa, subCNaC, S, M, Mc, N, ph, ph_arrivals)

    # Create the ODE matrices
    Φ, A, B = getFluidMatrices(ph, S)
    W = (Φ + B*P*A')

    # External arrivals
    u = vcat((p -> p.t).(ph_arrivals)...)
    λ = (x -> isempty(x) ? zeros(sum(C), 1) : Float64.(x))((E .* u)')

    # Reduced system for CN subset
    Pb = zeros(size(P))
    Pb[subCNa, subCNa] = P[subCNa, subCNa]
    Wb = (Φ + B*Pb*A')

    # Create queuing network struct
    QN = QueueNetworkStruct(Φ, A, B, P, W, E, λ, Pb, Wb, X0)

    return QN, Params     
end

function getQM(logfile, steps)
    matfile = matopen(joinpath(datapath, logfile))

    tv = read(matfile, "t")[:]
    qv = read(matfile, "q")[:]

    qsim = Matrix{Matrix{Float64}}(undef, size(tv, 1), size(qv[1], 2))
    for i = 1:size(qsim, 1)
        for j = 1:size(qsim, 2)
            if typeof(tv[i]) == Float64
                qsim[i, j] = [0.0 0.0]
            else
                qsim[i, j] = hcat(tv[i], qv[i][:, j])
            end
        end
    end

    qsteps = Matrix{Vector{Float64}}(undef, size(qsim))
    for (i, q) in enumerate(qsim)
        if size(q, 1) < 2
            qsteps[i] = q[1, 2] .* zeros(size(steps))
        else
            qsteps[i] = getQueueLengthAt(q, steps)
        end
    end
    #@assert all(all.([sum(qsteps[i, :]) .== sum(qsteps[1, :]) 
    #    for i = 1:size(qsteps, 1)]))

    return [mean(hcat(qsteps[:, i]...), dims=2)[:] for i=1:size(qsteps, 2)]
end

function extractData(suffix::String, Params::ParameterStruct; 
    T::Int64=0, printout=true, warn=true, get_alljobs=true)
    
    if printout; println("Extracting data for suffix $suffix"); end

    reqs = Array{Array{Float64, 2}, 1}(undef, sum(Params.C))
    alljobs = Dict{Int64, Array{Any, 2}}()
    jobs = Dict{String, Dict{Int64,Array{Any,2}}}()

    for (i, queue) in enumerate(Params.Queues)
        if Params.Disc[i] != "EXT"
            df_a = CSV.read(joinpath(datapath, "$(queue)_$suffix-Arv.csv"), DataFrame, silencewarnings=!warn)
            df_d = CSV.read(joinpath(datapath, "$(queue)_$suffix-Dep.csv"), DataFrame, silencewarnings=!warn)
            testData(df_a, df_d)

            # Reduce data    
            if T > 0
                df_a = df_a[df_a.TIMESTAMP .<= T, :]
                df_d = df_d[df_d.TIMESTAMP .<= T, :]
            end
            
            if printout; println("\tExtracting $queue"); end

            jobs[queue] = splitOnID(df_a, df_d)
        end
    end

    reqs_local = getRequestData(jobs, Params)

    n = 1
    idx = 1
    for (i, queue) in enumerate(Params.Queues)

        if Params.Disc[i] == "EXT"
            continue
        end

        for j in 1:length(Params.C_set[n])
            @assert all(diff(reqs_local[queue][Params.C_set[n][j]][:, 1]) .>= 0)
            @assert all(reqs_local[queue][Params.C_set[n][j]][:, 1] .<= 
                reqs_local[queue][Params.C_set[n][j]][:, 2])

            reqs[idx] =  reqs_local[queue][Params.C_set[n][j]]
            idx += 1
        end
        n += 1
    end

    if get_alljobs
        alljobs = getRequestIDData(jobs, Params)
    else
        alljobs = []
    end

    return getindex.(reqs, :, 1), getindex.(reqs, :, 2), alljobs
end

function testData(df_a::DataFrame, df_d::DataFrame)
    function testPerFrame(df::DataFrame)
        @assert !any(ismissing.(df.TIMESTAMP))
        @assert !any(ismissing.(df.JOB_ID))
        @assert !any(ismissing.(df.CLASS_ID))
        @assert all(diff(df.TIMESTAMP) .>= 0)
    end
    testPerFrame(df_a)
    testPerFrame(df_d)
end

function stackVisits(ta::DataFrame, td::DataFrame)

    visits = Array{Any,2}(undef,0,3)

    if size(ta, 1) == 0 && size(td, 1) == 0
        return visits
    end

    # Add first visit if queue starts non-empty
    i = findfirst(x -> x > ta.TIMESTAMP[1], td.TIMESTAMP)
    if !isnothing(i)
        @assert i <= 2 "job ID can max have 1 initial position!"
        if i-1 == 1
            visits = vcat(visits, [0.0 td.TIMESTAMP[1] td.CLASS_ID[1]])
        end
    else
        i = 1
    end

    # Add remaining matching between arrivals/departures
    Nd = size(td, 1) - size(visits,1)
    visits = vcat(visits, hcat(ta.TIMESTAMP[1:Nd], td.TIMESTAMP[i:end], td.CLASS_ID[i:end]))

    # Add remaining arrivals
    if size(ta, 1) - 2 <= Nd <= size(ta, 1) - 1
        visits = vcat(visits, hcat(ta.TIMESTAMP[Nd+1], 
                                [Inf], 
                                ta.CLASS_ID[Nd+1]))
    elseif  Nd < size(ta, 1) - 2
        error("Too many remaining arrivals per ID: $(Nd - size(ta, 1))")
    end

    # Assert that td larger than ta and each visit later than the previous, 
    #  and that their classes match
    @assert all(ta.CLASS_ID[1:Nd] .== td.CLASS_ID[i:end])
    @assert all(visits[:, 2] .>= visits[:, 1])
    @assert all(visits[2:end, 1] .>= visits[1:end-1, 2])

    return visits
end

function findIDmatch(df)
    m = maximum(df.JOB_ID)

    id_idx = [[] for k = 1:m+1]

    for k = 1:size(df, 1)
        append!(id_idx[df.JOB_ID[k]+1], k)
    end
    return id_idx
end

function splitOnID(df_a::DataFrame, df_d::DataFrame)
    IDs = union(unique(df_a.JOB_ID), unique(df_d.JOB_ID))
    jobs = Dict{Int64,Array{Any,2}}()
    id_idx_ta = findIDmatch(df_a)
    id_idx_td = findIDmatch(df_d)
    for id in IDs
        ta = df_a[id_idx_ta[id+1], :]
        td = id+1 > length(id_idx_td) ? df_d[[], :] : df_d[id_idx_td[id+1], :] 
        jobs[id] = stackVisits(ta, td)
    end
    return jobs
end

function getRequestData(jobs::Dict{String, Dict{Int64,Array{Any,2}}}, Params::ParameterStruct)
    
    reqsData = Dict{String, Dict{String, Array{Any, 2}}}()

    for (i, queue) in enumerate(Params.Queues)

        if Params.Disc[i] == "EXT"
            continue
        end

        ids = keys(jobs[queue])
 
        reqsData[queue] = Dict{String, Array{Any,2}}()
        for class in Params.Classes

            class_job_array = Array{Array{Float64, 2}, 1}(undef, length(ids))
            for (k, id) in enumerate(ids)
                jnk = jobs[queue][id]
                class_job_array[k] = jnk[jnk[:, 3] .== class, 1:2]
            end

            reqs = Array{Float64, 2}(undef, sum(size.(class_job_array, 1)), 2)
            reqs .= vcat(class_job_array...)
            idx = sortperm(reqs[:, 1])
            reqs = reqs[idx, :]

            @assert all(diff(reqs[:, 1]) .>= 0)
            @assert all(reqs[:, 2] .>= reqs[:, 1])
            @assert sum(.!isfinite.(reqs[:, 2])) <= length(ids)
            reqsData[queue][class] = reqs
        end

    end

    return reqsData
end

function getRequestIDData(jobs::Dict{String, Dict{Int64,Array{Any,2}}}, Params::ParameterStruct)
    alljobs = Dict{Int64, Array{Any,2}}()
    m = maximum(maximum.(keys.(values(jobs))))
    for id = 0:m
     
        n = size.((queue -> get(jobs[queue], id, [])).(Params.Queues[Params.Disc .!= "EXT"]), 1)
        alljobs[id] = Array{Any, 2}(undef, sum(n), 4)
        for (k, queue) in enumerate(Params.Queues[Params.Disc .!= "EXT"])
            if id in keys(jobs[queue])
                alljobs[id][(sum(n[1:k-1])+1):sum(n[1:k]), :] .= 
                    hcat(jobs[queue][id], repeat([queue], n[k], 1))
            end
        end 
        alljobs[id] = sortslices(alljobs[id], dims=1)


        @assert all(alljobs[id][1:end-1, 2] .== alljobs[id][2:end, 1])
        @assert all(diff(alljobs[id][:, 1]) .> 0)
    end
    return alljobs
end

function getQueueLengths(df_a::DataFrame, df_d::DataFrame, Q0::Dict{String,Int64})
    classes = union(unique(df_a.CLASS_ID), unique(df_d.CLASS_ID))
    q = Dict{String,Array{Float64,2}}()
    for class in classes
        ta = df_a[df_a.CLASS_ID .== class, :].TIMESTAMP
        td = df_d[df_d.CLASS_ID .== class, :].TIMESTAMP

        q[class] = zeros(length(ta) + length(td) + 1, 2)
        q[class][1, 2] = Q0[class]
        
        for k = 2:size(q[class], 1)
            # Findmin and pop from ta, td to increase/decrease queue length
            _, i = findmin([isempty(ta) ? Inf : ta[1], 
                            isempty(td) ? Inf : td[1]])
        
            if i == 2
                q[class][k, 1] = popfirst!(td)
                q[class][k, 2] = q[class][k-1, 2] - 1
            else
                q[class][k, 1] = popfirst!(ta)
                q[class][k, 2] = q[class][k-1, 2] + 1
            end
        end

        @assert all(q[class][:, 1] .>= 0)
        @assert all(q[class][:, 2] .>= 0)
        @assert all(-1 .<= diff(q[class][:, 2]) .<= 1)
        @assert all(diff(q[class][:, 1]) .>= 0)
    end
    
    if length(Classes) > 1
        q["All"] = addQueueLengths(q)
    end

    return q
end

function addQueueLengths(qDict::Dict{String,Array{Float64,2}})
    qv = collect(values(qDict))
    q = zeros(sum(size.(qv, 1)) - length(qv) + 1, 2)
    q[1, 2] = sum(getindex.(qv, 1, 2))
    idxs = ones(Int64, length(qv))
    for k = 2:size(q, 1)
        v = map((qc, idx) -> size(qc, 1) < idx+1 ?
                Inf : qc[idx+1, 1], qv, idxs)
        t, i = findmin(v)
        idxs[i] += 1
        q[k, 1] = t
        q[k, 2] = sum(getindex.(qv, idxs, 2)) 
    end 
    
    @assert all(q[:, 1] .>= 0)
    @assert all(q[:, 2] .>= 0)
    @assert all(-1 .<= diff(q[:, 2]) .<= 1)
    @assert all(diff(q[:, 1]) .>= 0)

    return q
end

function Q0ToDict(Q0::Array{Int64,2}, Queues::Array{String,1},
    Classes::Array{String,1})

    D = Dict{String,Dict{String, Int64}}()
    for (i, queue) in enumerate(Queues)
        D[queue] = Dict{String,Int64}()
        for (j, class) in enumerate(Classes)
            D[queue][class] = Q0[i,j]
        end
    end
    return D
end

function plotQueues(queues::Dict{String, Dict{String,Array{Float64,2}}}; 
      t_int=[0, 0], fignr=1)

    function plotQueue(q, class, t_int)

        if t_int[1] != t_int[2]
            t0 = max(1, findfirst(x -> x > t_int[1], q[:, 1]) - 1)
            tf = findfirst(x -> x > t_int[2], q[:, 1])
            tf = isnothing(tf) ? size(q, 1) : tf
        else
            t0 = 1
            tf = size(q, 1)
        end
        println(sum(q))
        println(t0)
        println(tf)
        step(q[t0:tf, 1], q[t0:tf, 2], where="post", label=class)
    end

    qkeys = sort(collect(keys(queues)))

    figure(fignr)
    clf()
    for (i, key) in enumerate(qkeys)
        subplot(length(qkeys), 1, i)

        classes = sort(collect(keys(queues[key])))
        if "All" in classes
            plotQueue(queues[key]["All"], "All", t_int)
        end
        for class in filter(s -> s != "All", classes)
            println(class)
            plotQueue(queues[key][class], class, t_int)
        end
        if t_int[1] != t_int[2]
            xlim(t_int)
        end
        title("$key, queue length")
        legend()
    end
end

# If Open, requires a single source/sink at the start/end
function getPHFromFile(csvfile)

    function getPH(X)
        v = Float64.(filter(x -> !ismissing(x), X))
        m = Int64(length(v)/2)

        Tv = v[1:m]
        tv = v[m+1:end]

        a = zeros(m)
        T = zeros(m, m)
        t = zeros(m)

        a[1] = 1.0
        for i = 1:m
            T[i, i] = -Tv[i]
            t[i] = Tv[i]*tv[i]

            if i < m; T[i, i+1] = Tv[i]*(1 - tv[i]); end
        end

        return EMpht.PhaseType(a, T, t, m)
    end

    ph_vals = CSV.read(csvfile, DataFrame, header=false)
    ph_source =  Array{EMpht.PhaseType, 1}(undef, sum(Disc[1] .== "EXT"))
    ph = Array{EMpht.PhaseType, 1}(undef, sum(C))

    if Disc[1] .== "EXT"
        ph_source[1] = getPH(Array(ph_vals[1, :]))
        for k = 1:length(ph)
            ph[k] = getPH(Array(ph_vals[k+1, :]))
        end
    else
        for k = 1:length(ph)
            ph[k] = getPH(Array(ph_vals[k, :]))
        end
    end

    return ph, ph_source
end