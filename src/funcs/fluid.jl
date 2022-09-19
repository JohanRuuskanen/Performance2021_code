

QNorm1(x::Array{Float64,1}, Params::ParameterStruct) = 
    ([sum(x[Params.M .== i]) for i = 1:Params.Q])[Params.M]

function g(x, c)
    return x > 0 ? min(c, x) / x : 1.0
end

function g_smooth(x, c, p)
    return 1/((1 + (x/c)^p)^(1/p))
end


function getFluidResults(p_opt::Array{Float64, 1}, QN::QueueNetworkStruct, Params::ParameterStruct;
        p_smooth=[])

    if isempty(p_smooth)
        p_smooth = p_opt
    end

    # Create the ODE problems
    alg = lsoda()
    x0 = QN.A*QN.X0
    
    prob_min = ODEProblem(dx_fluid!, x0, tspan, (QN, Params))
    prob_smooth = ODEProblem(dx_fluid_smooth!, x0, tspan, (QN, Params, p_smooth))
    prob_opt = ODEProblem(dx_fluid_smooth!, x0, tspan, (QN, Params, p_opt))

    # Extract the fluid solutions
    sol_min = solve(prob_min, alg, reltol=1e-8, abstol=1e-8)
    sol_smooth = solve(prob_smooth, alg, reltol=1e-8, abstol=1e-8)
    sol_opt = solve(prob_opt, alg, reltol=1e-8, abstol=1e-8)

    qf_min = [sum(hcat(sol_min.u...)[Params.N .== i, :], dims=1)[:] for i in 1:sum(Params.C)]
    qf_smooth = [sum(hcat(sol_smooth.u...)[Params.N .== i, :], dims=1)[:] for i in 1:sum(Params.C)]
    qf_opt = [sum(hcat(sol_opt.u...)[Params.N .== i, :], dims=1)[:] for i in 1:sum(Params.C)]

    # Extract the fixed points per class
    qm_min = getindex.(qf_min, length.(qf_min))
    qm_smooth = getindex.(qf_smooth, length.(qf_smooth))
    qm_opt = getindex.(qf_opt, length.(qf_opt))

    # Extract the estimated processor share for both min and smooth
    ps_min = g.(QNorm1(sol_min.u[end], Params), Params.K[Params.M])
    ps_smooth = g_smooth.(QNorm1(sol_smooth.u[end], Params), Params.K[Params.M], p_smooth[Params.M])
    ps_opt = g_smooth.(QNorm1(sol_opt.u[end], Params), Params.K[Params.M], p_opt[Params.M])

    # Extract the functions for RT CDF approximation
    rt_cdf_min = Array{Function, 1}(undef, sum(Params.C))
    rt_cdf_smooth = Array{Function, 1}(undef, sum(Params.C))
    rt_cdf_opt = Array{Function, 1}(undef, sum(Params.C))

    for i = 1:sum(Params.C)
        j = findfirst(x -> i <= x, cumsum(Params.C))

        rt_cdf_min[i] = (t -> tr_cdf(t, Params.ph[i].T', sol_min.u[end][Params.M .== j], 
            Params.K[j], Params.ph[i].π))
        rt_cdf_smooth[i] = (t -> tr_cdf_smooth(t, Params.ph[i].T', sol_smooth.u[end][Params.M .== j], 
            Params.K[j], p_smooth[j], Params.ph[i].π))
        rt_cdf_opt[i] = (t -> tr_cdf_smooth(t, Params.ph[i].T', sol_opt.u[end][Params.M .== j], 
            Params.K[j], p_opt[j], Params.ph[i].π))
    end 

    # Extract the functions for RT CDF approximation over CN
    β_min = getBetaInData(sol_min.u[end], ps_min, QN, Params)
    β_smooth = getBetaInData(sol_smooth.u[end], ps_smooth, QN, Params)
    β_opt = getBetaInData(sol_opt.u[end], ps_opt, QN, Params)

    rt_cdf_CN_min(t) = tr_cdf_subCN(t, β_min, ps_min, QN)
    rt_cdf_CN_smooth(t) = tr_cdf_subCN(t, β_smooth, ps_smooth, QN)
    rt_cdf_CN_opt(t) = tr_cdf_subCN(t, β_opt, ps_opt, QN)

    return FluidResultStruct(sol_min, qf_min, qm_min, ps_min, rt_cdf_min, rt_cdf_CN_min,
        sol_smooth, qf_smooth, qm_smooth, ps_smooth, p_smooth, rt_cdf_smooth, rt_cdf_CN_smooth,
        sol_opt, qf_opt, qm_opt, ps_opt, p_opt, rt_cdf_opt, rt_cdf_CN_opt)
end

function fluidSol(Trace, QN, Params, ts, steps)
    alg = lsoda()
    x0 = QN.A*QN.X0

    prob_min = ODEProblem(dx_fluid!, x0, float(ts), (QN, Params))
    prob_smooth = ODEProblem(dx_fluid_smooth!, x0, float(ts), (QN, Params, p_smooth))
    prob_opt = ODEProblem(dx_fluid_smooth!, x0, float(ts), (QN, Params, Trace.p_opt))

    # Extract the fluid solutions
    sol_min = solve(prob_min, alg, saveat=steps, reltol=1e-8, abstol=1e-8)
    sol_smooth = solve(prob_smooth, saveat=steps, alg, reltol=1e-8, abstol=1e-8)
    sol_opt = solve(prob_opt, alg, saveat=steps, reltol=1e-8, abstol=1e-8)

    qf_min = [sum(hcat(sol_min.u...)[Params.N .== i, :], dims=1)[:] for i in 1:sum(Params.C)]
    qf_smooth = [sum(hcat(sol_smooth.u...)[Params.N .== i, :], dims=1)[:] for i in 1:sum(Params.C)]
    qf_opt = [sum(hcat(sol_opt.u...)[Params.N .== i, :], dims=1)[:] for i in 1:sum(Params.C)]

    return qf_min, qf_smooth, qf_opt

end


function getFluidMatrices(ph::Array{EMpht.PhaseType, 1}, S::Array{Int64,1})
    Φ = blockdiag(sparse.((x->x.T).(ph))...)
    A = blockdiag(sparse.(reshape.((x->x.π).(ph), S, 1))...)
    B = blockdiag(sparse.(reshape.((x->x.t).(ph), S, 1))...)
    return Φ, A, B
end

# p[1]: QueueNetworkStruct, p[2]: ParameterStruct
function dx_fluid!(dx, x, p, t)
    dx .= p[1].W' * (g.(QNorm1(x, p[2]), p[2].K[p[2].M]) .* x) + sum(p[1].A*p[1].λ, dims=2)[:]
end

# p[1]: QueueNetworkStruct, p[2]: ParameterStruct, p[3]: p_smooth values
function dx_fluid_smooth!(dx, x, p, t)
    dx .= p[1].W' * (g_smooth.(QNorm1(x, p[2]), p[2].K[p[2].M], p[3][p[2].M]) .* x) + 
        sum(p[1].A*p[1].λ, dims=2)[:]
end

tr_cdf(t, Φ, x0, C, α) = 1 .- sum(exp((min(C, sum(x0))/sum(x0))*Φ*t)*α)
tr_cdf_smooth(t, Φ, x0, C, p, α) = 1 .- sum(exp(1/norm([1, sum(x0)/C], p)*Φ*t)*α)
tr_cdf_qdist(t, cf::Function, qdist) = sum([cf(t, k)*qdist[k] for k = 1:length(qdist)])

function tr_cdf_subCN(t::Float64, β::Array{Float64, 1}, ps::Array{Float64, 1}, QN::QueueNetworkStruct)
    return (1 .- (QN.A*β)' * exp(Diagonal(ps) * QN.Wb * t) * ones(size(QN.Φ, 1), 1))[1]
end

function tr_cdf_subCN_dist(t::Float64, β::Array{Float64, 1}, qds::Array{Float64, 2}, 
        QN::QueueNetworkStruct, Params::ParameterStruct)
    p = 0
    N = size(qds, 1)
    for k in 1:N
        p += (1/N * (QN.A*β)' * exp(Diagonal(g.((qds[k, :] .+ 1)[Params.M], Params.K[Params.M])) * 
            QN.Wb * t) * ones(size(QN.Φ, 1), 1))[1]
    end
    return 1 - p
end

function getBetaInData(xs::Array{Float64, 1}, ps::Array{Float64, 1}, 
        QN::QueueNetworkStruct, Params::ParameterStruct)
    Eb = QN.P[Params.subCNaC, Params.subCNa]
    Fb = getOutFlow(xs, ps, QN)
    β_flow = zeros(sum(Params.C))
    β_flow[Params.subCNa] .= normalize(Eb' * Fb[Params.subCNaC], 1)
    return β_flow
end

function tr_cdf_data(t, trs, pvec)
    if t < 0; return 0.0; end
    idx = findfirst(x -> x > t, trs)
    return isnothing(idx) ? 1.0 : pvec[idx]
end

function findQuantile(f, α; x0=1.0, order=0)
    return fzero((t -> f(t) - α), x0, order=order)
end

function findQuantileData(t, α)
    ts = sort(t)
    pvec = (1:length(ts)) ./ length(ts)
    k = 1
    while pvec[k] < α
        if k == length(ts)
            break
        end
        k += 1
    end
    return ts[k]
end