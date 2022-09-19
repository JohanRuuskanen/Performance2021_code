
tryparse(::Type{Int64}, a::Int64) = a

mse(e) = mean(e.^2)
mae(e) = maximum(abs.(e))


removeNothing(x::Array{T, 1}) where {T <: Any} = x[.!isnothing.(x)]
function removeNothing(x::Array{Union{Nothing, T}, 1}) where {T <: Any} 
    return Array{T, 1}(x[.!isnothing.(x)])
end

nothing2Missing(x) = isnothing(x) ? missing : x
function nothing2Missing(x::AbstractArray{T, 1}, i::Union{Nothing, Integer}) where {T <: Any}
    return isnothing(i) ? missing : x[i]
end
function nothing2Missing(f::Function, x::AbstractArray{T, 1}, 
        i::Union{Nothing, Integer}) where {T <: Any}
    return isnothing(i) ? missing : f(x[i])
end

function getTraceData(suffix::String, QN::QueueNetworkStruct, Params::ParameterStruct; 
        printout=true, warn=true, T=0, get_for_CN=true)
    # Extract all arrivals/departures
    ta, td, alljobs = extractData(suffix, Params, 
        T=T, printout=printout, warn=warn, get_alljobs=get_for_CN)

    @assert all([all(ta[k] .>= 0) for k = 1:sum(Params.C)])
    @assert all([all(td[k] .>= 0) for k = 1:sum(Params.C)])
    @assert all([all(td[k] .>= ta[k]) for k = 1:sum(Params.C)])

    # Calculate queue lengths per class and per queue
    q = getQueueLengths.(ta, td)
    q_node = [addQueueLengths(q[Params.C_idx[k]]) for k in 1:Params.Q]

    # Get queue length distributio
    qd = getqdist.(q_node)
    @assert all(isapprox.(sum.(qd), 1))

    # Extract response times per class
    tw = filter.(x -> isfinite(x), td - ta)
   
    @assert all(all.((t -> t .>= 0).(tw)))

    # Get the steady state queue lengths per class and per queue
    qm = [mean(q[k][1:end-1, 2], weights(diff(q[k][:, 1]))) for k = 1:sum(Params.C)]
    qm_node = [mean(q_node[k][1:end-1, 2], weights(diff(q_node[k][:, 1]))) for k = 1:Params.Q]

    # Extract optimal p-norm parameter and processor share per queue
    p_opt = getOptimPNorm(q_node, Params)
    ps = getAvgPS.(q_node, Params.K)

    # Extract values for CN
    if get_for_CN
        taCN, tdCN = getDataOverN(alljobs, Params.subCN)

        @assert all(taCN .>= 0)
        @assert all(tdCN .>= 0)
        @assert all(tdCN .>= taCN)

        # Get queue lengths over CN subset
        qCN = getQueueLengths(taCN, tdCN)
        qa = addQueueLengths(q[Params.subCNa])
        @assert sum(qa - qCN) == 0.0

        # Extract response times over CN subset
        twCN = filter(x -> isfinite(x), tdCN - taCN)
        @assert all(twCN .>= 0)
    else
        taCN = Array{Float64, 1}(undef, 0)
        tdCN = Array{Float64, 1}(undef, 0)
        twCN = Array{Float64, 1}(undef, 0)
        qCN = Array{Float64, 2}(undef, 0, 2)
    end

    return TraceStruct(ta, td, taCN, tdCN, q, q_node, qd, qCN, tw, twCN, qm, qm_node, p_opt, ps)
end

function getQueueLengths(ta_orig::Array{Float64,1}, td_orig::Array{Float64,1})

    @assert size(ta_orig) == size(td_orig)
    @assert all(ta_orig .<= td_orig)

    ta = sort(ta_orig)
    td = sort(td_orig)

    m = length(td)
    q = zeros(2*m+1, 2)

    ka = 1
    kd = 1
    for k = 2:size(q, 1)
        # Findmin and pop from ta, td to increase/decrease queue length
        v, i = findmin([ka > m ? Inf : ta[ka], 
                        kd > m ? Inf : td[kd]])

        if !isfinite(v)
            q = q[1:k-1, :]
            break
        end

        if i == 2
            q[k, 1] = v
            q[k, 2] = q[k-1, 2] - 1
            kd += 1
        else
            q[k, 1] = v
            q[k, 2] = q[k-1, 2] + 1
            ka += 1
        end

    end

    @assert all(q[:, 1] .>= 0)
    @assert all(q[:, 2] .>= 0)
    @assert all(-1 .<= diff(q[:, 2]) .<= 1)
    @assert all(diff(q[:, 1]) .>= 0)

    return q
end

function getQueueLengths(ta::Array{DateTime,1}, td::Array{DateTime,1})
    q = getQueueLengths(datetime2unix.(ta), datetime2unix.(td))
    return hcat(unix2datetime.(q[:,1]), q[:, 2])
end

function movingAvg(x::Array{T, 1}, w::AbstractWeights; 
        window::Int64=100, α=0.95) where T <: Number

    idxf = findfirst(x -> x > 0, w)
    w[1:idxf] .= 1e-6
    
    @assert length(x) == length(w)
    xm = zeros(size(x))
    xa = zeros(size(x))
    for k = 1:length(x)
        if k <= window
            xm[k] = mean(view(x, 1:k), w[1:k])
            #xa[k] = quantile(view(x, 1:k), w[1:k], α)
        else
            xm[k] = mean(view(x, k-window:k), w[k-window:k])
            #xa[k] = quantile(view(x, k-window:k), w[k-window:k], α)
        end
    end
    return xm, xa
end

function movingAvg(x::Array{T, 1}; window::Int64=100, α=0.95) where T <: Number
    xm = zeros(size(x))
    xa = zeros(size(x))
    for k = 1:length(x)
        if k <= window
            xm[k] = mean(view(x, 1:k))
            xa[k] = quantile(view(x, 1:k), α)
        else
            xm[k] = mean(view(x, k-window:k))
            xa[k] = quantile(view(x, k-window:k), α)
        end
    end
    return xm, xa
end

function addQueueLengths(q_v::Array{Array{Float64, 2}, 1})
    q_tot = zeros(sum(size.(q_v, 1)) - length(q_v) + 1, 2)
    idxs = ones(Int64, length(q_v))
    for k = 2:size(q_tot,1)
        vals = map((q, idx) -> size(q, 1) < idx+1 ?
                Inf : q[idx+1, 1], q_v, idxs)

        v, i = findmin(vals)
        idxs[i] += 1
        q_tot[k, 1] = v
        q_tot[k, 2] = sum(getindex.(q_v, idxs, 2))
    end

    idxs = ones(Bool, size(q_tot, 1))
    for k = 2:size(q_tot, 1)-1
        if idxs[k] == 0
            continue
        end

        if q_tot[k, 1] == q_tot[k+1, 1] && q_tot[k-1, 2] == q_tot[k+1, 2]
            idxs[k:k+1] .= 0
        end
    end

    q_tot = q_tot[findall(idxs), :]

    @assert all(q_tot[:, 1] .>= 0)
    @assert all(q_tot[:, 2] .>= 0)
    @assert all(-1 .<= diff(q_tot[:, 2]) .<= 1)
    @assert all(diff(q_tot[:, 1]) .>= 0)

    return q_tot
end

function getAvgQueueLengths(ta::Array{Float64, 1}, td::Array{Float64, 1})
    q = getQueueLengths(ta, td)
    qa = zeros(size(ta))

    ka = 1
    for i in 1:length(qa)
        while q[ka, 1] != ta[i]
            ka += 1
            if q[ka, 1] > ta[i]
                error("ta[i] could not be found in q!")
            end
        end

        kd = ka
        while q[kd, 1] != td[i]
            kd += 1
            if kd > size(q, 1)
                error("td[i] could not be found in q!")
            end
            qa[i] += (q[kd, 1] - q[kd-1, 1])*q[kd-1, 2]
        end
        qa[i] /= td[i] - ta[i]
    end

    return qa
end

function getAvgQueueLengthsOverSim(t::AbstractArray{Float64, 1})

    q_avg = Array{Array{Float64, 2}, 1}(undef, sum(C))
    for k = 1:sum(C)
        q_avg[k] = hcat(t, zeros(length(t)))
    end

    for sim = 1:MCsims
        ta, td, alljobs = extractData(sim, T=round(Int64, maximum(t)), warn=false, printout=false)
        q = getQueueLengths.(ta, td)
        for k = 1:sum(C)
            q_avg[k][:, 2] += getQueueLengthAt(q[k], t)
        end
    end

    for k = 1:sum(C)
        q_avg[k] = q_avg[k] ./ repeat([1.0 MCsims], length(t), 1)
    end

    return q_avg
end

function getQueueLengthAt(q::Array{Float64, 2}, t::AbstractArray{T, 1}) where T <: Number

    function re_sort(q, idx)
        q_t = zeros(size(q))
        q_t[idx_s] = q
        return q_t
    end

    idx_s = sortperm(t)
    t_sort = t[idx_s]
    q_t_sort = zeros(length(t))

    k = 2
    for i = 1:length(t)
        while t_sort[i] >= q[k, 1]
            k += 1
            if k > size(q, 1)
                q_t_sort[i:end] .= q[end, 2]
                return re_sort(q_t_sort, idx_s)
            end
        end
        q_t_sort[i] = q[k-1, 2]
    end

    q_d = re_sort(q_t_sort, idx_s)
    @assert all(q_d .>= 0)
    
    return q_d
end

function getUtil(q::Array{Float64, 2}, k::T) where T <: Number

    if k == Inf
        return 0.0
    else
        k == round(Int64, k)
    end
    
    t0 = 0
    for i = 2:size(q, 1)
        for j = 0:k-1
            if q[i-1, 2] == j
                t0 += (k-j)*(q[i, 1] - q[i-1, 1])
            end
        end
    end

    return 1 - (t0 / (k*q[end, 1]))
end

function getUtil(ta::Array{Float64, 1}, td::Array{Float64, 1}, k::T; dt=1, tw=10) where T <: Number

    if k == Inf
        return 0.0
    else
        k == round(Int64, k)
    end
    
    q = getQueueLengths(ta, td)

    u = q[2,1]:dt:q[end,1]
    u = hcat(u, zeros(size(u)))
    
    j = 1
    for i in 2:size(u, 1)
        while q[j, 1] < u[i, 1]
            j += 1
            if j > size(q, 1)
                println("$(u[i, 1]), $(q[j-1, 1])")
                error("j larger than size of q!")
            end
        end

        # queue length at u[i, 1] is q[j-1, 2]
        u[i, 2] += min(k, q[j-1, 2])*(u[i, 1] - q[j-1, 1])

        t = u[i, 1] - q[j-1, 1]
        l = 1

        while t + q[j-l, 1] - q[j-l-1, 1] < tw
            u[i, 2] += min(k, q[j-l-1, 2])*(q[j-l, 1] - q[j-l-1, 1])
            t += q[j-l, 1] - q[j-l-1, 1]
            l += 1
        end

        if j > l+1
            u[i, 2] += min(k, q[j-l-1, 2])*(q[j-l, 1] - (u[i, 1] - tw))
        end
        t += q[j-l, 1] - (u[i, 1] - tw)

        u[i, 2] /= t

        @assert t == tw

    end

    return u
end

function getUsedNodes(P::Array{Float64, 2}, Q0::Array{Int64, 2})
    idx_sources = (sum(P, dims=1)[:] .== 0) .& (sum(P, dims=2)[:] .>0)
    idx_sinks = (sum(P, dims=1)[:] .> 0) .& (sum(P, dims=2)[:] .== 0)
    idx_trans = (sum(P, dims=1)[:] .> 0) .& (sum(P, dims=2)[:] .> 0)
    idx_start = Q0'[:] .> 0
    return idx_trans .| idx_start, idx_sources, idx_sinks
end

function classQueueP(Pall::Array{Float64, 2}, Queues, Classes)
    Q = length(Queues)
    C = length(Classes)
    
    P = zeros(Q, Q, C, C)
    for k = 1:C
        for i = 1:C
            mv = ((k-1)*Q + 1):k*Q
            nv = ((i-1)*Q + 1):i*Q
            P[:, :, k, i] = Pall[mv, nv]
        end
    end
    return P
end

function P_NCP_splat(P::Array{Float64, 4})
    m, _, n, _ = size(P)
    @assert size(P) == (m, m, n, n)

    P_new = zeros(m*n, m*n)
    for i = 1:m
        for j = 1:m
            idx_r = (i-1)*n+1:(i-1)*n+n
            idx_c = (j-1)*n+1:(j-1)*n+n
            P_new[idx_r, idx_c] = P[i, j, :, :]
        end
    end
    return P_new
end

function getAvgPS(q::Array{Float64, 2}, k)
    idx_gk = findall(x -> x > k, q[1:end-1, 2])
    idx_0 =  findall(x -> x == 0, q[1:end-1, 2])

    T = q[end, 1] - q[1, 1]
    ts = diff(q[:, 1])

    return (T - sum(ts[idx_0]) - sum(ts[idx_gk]))/T + sum(k ./ q[idx_gk, 2] .* ts[idx_gk])/T
end

function getSmooth(q::Array{Float64, 2}, k)
    idx_k =  findall(x -> 1 <= x <= k, q[1:end-1, 2])
    idx_gk =  findall(x -> x > k, q[1:end-1, 2])

    T = q[end, 1] - q[1, 1]
    ts = diff(q[:, 1])

    return sum(ts[idx_k].*q[idx_k, 2])/T + k*sum(ts[idx_gk])/T
end

function getqdist(q::Array{Float64, 2})
    m = Int64(maximum(q[:, 2])) + 1

    if m == 1
        return [1.0]
    end

    qd = zeros(m)

    for k = 1:size(q, 1)-1
        qd[Int64(q[k, 2])+1] += q[k+1, 1] - q[k, 1]
    end

    qd ./= q[end, 1] - q[1, 1]

    return qd
end

function getDataOverN(alljobs::Dict{Int64, Array{Any, 2}}, N::Array{String, 2})
    
    Nv = [N[k, :] for k=1:size(N, 1)]
    tdN = -1*ones(sum(size.(values(alljobs), 1)))
    taN = -1*ones(sum(size.(values(alljobs), 1)))
    c = 1

    for key in keys(alljobs)
        k = 1
        while k <= size(alljobs[key], 1)
            if alljobs[key][k, 3:4] in Nv
                j = k
                while alljobs[key][j, 3:4] in Nv
                    j += 1
                    if j > size(alljobs[key], 1)
                        break
                    end
                end

                tdN[c] = alljobs[key][j-1, 2]
                taN[c] = alljobs[key][k, 1]   
                c += 1
                
                k = j
            end
            k += 1
        end
    end

    idx = sortperm(taN[1:c-1])
    taN = taN[idx]
    tdN = tdN[idx]

    return taN, tdN
end

function getOptimPNorm(q_node::Array{Array{Float64, 2}, 1}, Params::ParameterStruct)
    p_opt = zeros(Params.Q)
    for k = 1:Params.Q
        qnm = mean(q_node[k][1:end-1, 2], weights(diff(q_node[k][:, 1])))
        if Params.Disc[Params.Disc .!= "EXT"][k] == "INF"
            p_opt[k] = Inf
        elseif Params.Disc[Params.Disc .!= "EXT"][k] == "PS"
            try
                p_opt[k] = fzero((p -> qnm/norm([Params.K[k], qnm], p) - getUtil(q_node[k], Int64(Params.K[k]))), 1.0)
            catch e
                println("Warning: fzero got error $e, setting p_opt 0")
                p_opt[k] = 0.0
            end
        end
    end
    return p_opt
end

function getOptimPNorm(q_node::Array{Array{Float64, 2}, 1}, ρ::Array{Float64, 1})
    p_opt = zeros(Q)
    for k = 1:Q
        qnm = mean(q_node[k][1:end-1, 2], weights(diff(q_node[k][:, 1])))
        if Disc[k] == "INF"
            p_opt[k] = Inf
        elseif Disc[k] == "PS"
            try
                p_opt[k] = fzero((p -> qnm/norm([K[k], qnm], p) - ρ[k]), 0.1)
            catch
                p_opt[k] = 0.0
            end
        end
    end
    return p_opt
end

function getPhasePos(S::Array{Int64,1}, Q::Int64, C::Array{Int64,1})

    # Which state belongs to which queue
    M = vcat([i*ones(Int64, sum((S[sum(C[1:i-1])+1:sum(C[1:i])]))) for i in 1:Q]...)

    # Which class belongs to which queue
    Mc = vcat([i*ones(Int64, C[i]) for i in 1:Q]...)

    # Which state belongs to which class
    N = vcat([i*ones(Int64, S[i]) for i in 1:sum(C)]...)

    return M, Mc, N
end

function getPhasesPos(ph::Array{EMpht.PhaseType, 1}, Q::Int64, C::Array{Int64,1})
    # Phases per class/queue
    S = (x -> x.p).(ph)

    @assert length(S) == sum(C)
    
    M, Mc, N = getPhasePos(S, Q, C)
    return S, M, Mc, N
end

function sampleQD(qd, samp)
    qds = zeros(samp, length(qd))
    for k = 1:length(qd)
        qds[:, k] = sample(0:length(qd[k])-1, weights(qd[k]), samp)
    end
    return qds
end

function getOutFlow(xs::Array{Float64, 1}, ps::Array{Float64, 1}, QN::QueueNetworkStruct)
    return QN.B' * (ps .* xs)
end

function getAvgRT(xs::Array{Float64, 1}, ps::Array{Float64, 1}, 
        QN::QueueNetworkStruct, Params::ParameterStruct)
    return [sum(xs[Params.N .== i]) for i = 1:sum(Params.C)] ./ getOutFlow(xs, ps, QN)
end

function getAvgRTCN(xs::Array{Float64, 1}, ps::Array{Float64, 1},  
        QN::QueueNetworkStruct, Params::ParameterStruct; Itrs=100)
    trm = getAvgRT(xs, ps, QN, Params)
    β = getBetaInData(xs, ps, QN, Params)
    return  β' * (sum([QN.Pb^k for k = 0:Itrs])) * trm
end

function joinDataOverNode(ta::Array{Array{Float64, 1}, 1}, td::Array{Array{Float64, 1}, 1})
    sort_node = Array{Array{Int64, 1}, 1}(undef, Q)
    ta_node = Array{Array{Float64, 1}, 1}(undef, Q)
    td_node = Array{Array{Float64, 1}, 1}(undef, Q)
    for i = 1:Q
        idx = findall(Mc .== i)
        ta_node[i] = vcat(ta[idx]...)
        td_node[i] = vcat(td[idx]...)

        idx_sort = sortperm(ta_node[i])
        sort_node[i] = idx_sort
        ta_node[i] = ta_node[i][idx_sort]
        td_node[i] = td_node[i][idx_sort]
    end
    return ta_node, td_node, sort_node
end

function getPHFromMatrix(ph_params, nbrC, hasSource)

    function getPH(X)
        m = length(X[1])

        if m == 1
            return EMpht.PhaseType([1.0], repeat([-X[1][1]], 1, 1), [X[1][1]], m)
        end

        a = zeros(m)
        T = zeros(m, m)
        t = zeros(m)

        a[1] = 1.0
        for i = 1:m
            T[i, i] = -X[1][i]
            t[i] = X[1][i]*X[2][i]

            if i < m; T[i, i+1] = X[1][i]*(1 - X[2][i]); end
        end

        return EMpht.PhaseType(a, T, t, m)
    end

    ph_arrivals =  Array{EMpht.PhaseType, 1}(undef, hasSource)
    ph = Array{EMpht.PhaseType, 1}(undef, nbrC)

    if hasSource
        ph_arrivals[1] = getPH(ph_params[1, :])
        for k = 1:length(ph)
            ph[k] = getPH(Array(ph_params[k+1, :]))
        end
    else
        for k = 1:length(ph)
            ph[k] = getPH(Array(ph_params[k, :]))
        end
    end

    return ph, ph_arrivals
end

mp(x, xh, f) = f(abs.(x - xh) ./ x)

function compareResults(Results::FluidResultStruct, Trace::TraceStruct; a=0.95)

    function tryquadk(f::Function, suffix::String)
        try
            return quadgk(f, 0, Inf)[1]
        catch
            println("Warning, quadfk failed on $suffix")
            return -1.0
        end
    end

    function tryfzero(f::Function, suffix::String)
        try
            return fzero(f, 1.0, order=0)
        catch
            println("Warning, fzero failed on $suffix")
            return -1.0
        end
    end

    q_meanp_min = mp(Trace.qm, Results.qm_min, mean) 
    q_meanp_smooth = mp(Trace.qm,Results.qm_smooth, mean)
    q_meanp_opt = mp(Trace.qm, Results.qm_opt, mean)

    q_maxp_min = mp(Trace.qm, Results.qm_min, maximum) 
    q_maxp_smooth = mp(Trace.qm,Results.qm_smooth, maximum)
    q_maxp_opt = mp(Trace.qm, Results.qm_opt, maximum)

    trm = zeros(length(Trace.qm), 4)
    tra = zeros(length(Trace.qm), 4)

    trmCN = zeros(4)
    traCN = zeros(4)

    for i = 1:length(Trace.qm)
        trm[i, 1] = mean(Trace.tw[i])
        trm[i, 2] = tryquadk(t -> 1 - Results.rt_cdf_min[i](t), "min")
        trm[i, 3] = tryquadk(t -> 1 - Results.rt_cdf_smooth[i](t), "min")
        trm[i, 4] = tryquadk(t -> 1 - Results.rt_cdf_opt[i](t), "min")

        tra[i, 1] = findQuantileData(Trace.tw[i], a)
        tra[i, 2] = tryfzero(t-> Results.rt_cdf_min[i](t) - a, "min")
        tra[i, 3] = tryfzero(t-> Results.rt_cdf_smooth[i](t) - a, "smooth")
        tra[i, 4] = tryfzero(t-> Results.rt_cdf_opt[i](t) - a, "opt")
    end

    trmCN[1] = mean(Trace.twCN)
    trmCN[2] = tryquadk(t -> 1 - Results.rt_cdf_CN_min(t), "min")
    trmCN[3] = tryquadk(t -> 1 - Results.rt_cdf_CN_smooth(t), "smooth")
    trmCN[4] = tryquadk(t -> 1 - Results.rt_cdf_CN_opt(t), "opt")

    traCN[1] = findQuantileData(Trace.twCN, a)
    traCN[2] = tryfzero(t-> Results.rt_cdf_CN_min(t) - a, "min")
    traCN[3] = tryfzero(t-> Results.rt_cdf_CN_smooth(t) - a, "smooth")
    traCN[4] = tryfzero(t-> Results.rt_cdf_CN_opt(t) - a, "opt")

    trm_meanp_min = mp(trm[:, 1], trm[:, 2], mean)
    trm_meanp_smooth = mp(trm[:, 1], trm[:, 3], mean)
    trm_meanp_opt = mp(trm[:, 1], trm[:, 4], mean)

    trm_maxp_min = mp(trm[:, 1], trm[:, 2], maximum)
    trm_maxp_smooth = mp(trm[:, 1], trm[:, 3], maximum)
    trm_maxp_opt = mp(trm[:, 1], trm[:, 4], maximum)

    tra_meanp_min = mp(tra[:, 1], tra[:, 2], mean)
    tra_meanp_smooth = mp(tra[:, 1], tra[:, 3], mean)
    tra_meanp_opt = mp(tra[:, 1], tra[:, 4], mean)

    tra_maxp_min = mp(tra[:, 1], tra[:, 2], maximum)
    tra_maxp_smooth = mp(tra[:, 1], tra[:, 3], maximum)
    tra_maxp_opt = mp(tra[:, 1], tra[:, 4], maximum)

    trm_CN_min = mp(trmCN[1], trmCN[2], x->x)
    trm_CN_smooth = mp(trmCN[1], trmCN[3], x->x)
    trm_CN_opt = mp(trmCN[1], trmCN[4], x->x)

    tra_CN_min = mp(traCN[1], traCN[2], x->x)
    tra_CN_smooth = mp(traCN[1], traCN[3], x->x)
    tra_CN_opt = mp(traCN[1], traCN[4], x->x)

    return ErrorStruct(
        q_meanp_min,
        q_meanp_smooth,
        q_meanp_opt,
        q_maxp_min,
        q_maxp_smooth,
        q_maxp_opt,
        trm_meanp_min,
        trm_meanp_smooth,
        trm_meanp_opt,
        trm_maxp_min,
        trm_maxp_smooth,
        trm_maxp_opt,
        tra_meanp_min,
        tra_meanp_smooth,
        tra_meanp_opt,
        tra_maxp_min,
        tra_maxp_smooth,
        tra_maxp_opt,
        trm_CN_min,
        trm_CN_smooth,
        trm_CN_opt,
        tra_CN_min,
        tra_CN_smooth,
        tra_CN_opt
    )

end

function printErrors(Errors::ErrorStruct)
    println("Steady-state errors in percent")
    println("Metric, standard, smoothed, optimal")
    println("q_mean: $(round.(Errors.q_meanp_min, digits=2)) $(round.(Errors.q_meanp_smooth, digits=2)) "*
        "$(round.(Errors.q_meanp_opt, digits=2))")
    println("q_max: $(round.(Errors.q_maxp_min, digits=2)) $(round.(Errors.q_maxp_smooth, digits=2)) "*
        "$(round.(Errors.q_maxp_opt, digits=2))")

    println("trm_mean: $(round.(Errors.trm_meanp_min, digits=2)) $(round.(Errors.trm_meanp_smooth, digits=2)) "*
        "$(round.(Errors.trm_meanp_opt, digits=2))")
    println("trm_max: $(round.(Errors.trm_maxp_min, digits=2)) $(round.(Errors.trm_maxp_smooth, digits=2)) "*
        "$(round.(Errors.trm_maxp_opt, digits=2))")

    println("tra_mean: $(round.(Errors.tra_meanp_min, digits=2)) $(round.(Errors.tra_meanp_smooth, digits=2)) "*
        "$(round.(Errors.tra_meanp_opt, digits=2))")
    println("tra_max: $(round.(Errors.tra_maxp_min, digits=2)) $(round.(Errors.tra_maxp_smooth, digits=2)) "*
        "$(round.(Errors.tra_maxp_opt, digits=2))")

    println("")
    println("trm_CN: $(round.(Errors.trm_CN_min, digits=2)) $(round.(Errors.trm_CN_smooth, digits=2)) "*
        "$(round.(Errors.trm_CN_opt, digits=2))")
    println("tra_CN: $(round.(Errors.tra_CN_min, digits=2)) $(round.(Errors.tra_CN_smooth, digits=2)) "*
        "$(round.(Errors.tra_CN_opt, digits=2))")

end