function tr_cdf_plot_data(x, tw::Array{Float64, 1}, α::Float64)
    trs = sort(tw)
    pvec = (1:length(tw)) ./ length(tw)
    trm = mean(tw)
    tra = findQuantileData(tw, α)

    y = (t->tr_cdf_data(t, trs, pvec)).(x)
    ym = NaN*ones(size(x))
    ya = NaN*ones(size(x))
    ym[findfirst(x -> x > trm, x)] = tr_cdf_data(trm, trs, pvec)
    ya[findfirst(x -> x > tra, x)] = tr_cdf_data(tra, trs, pvec)

    plot(x, y, "C0", label="data")
    plot(trm, tr_cdf_data(trm, trs, pvec), "C0o")
    plot(tra, tr_cdf_data(tra, trs, pvec), "C0^")

    return hcat(y, ym, ya)
end

function tr_cdf_plot(x, cf::Function, α::Float64, color::String, labelstr::String)
    y = cf.(x)
    ym = NaN*ones(size(x))
    ya = NaN*ones(size(x))

    plot(x, y, color, label=labelstr)
    try
        trm = quadgk((t -> 1 - cf(t)), 0, Inf)[1]
        tra = fzero((t -> cf(t) - α), 1.0, order=0)

        ym[findfirst(x -> x > trm, x)] = cf(trm)
        ya[findfirst(x -> x > tra, x)] = cf(tra)

        plot(trm, cf(trm), color*"o")
        plot(tra, cf(tra), color*"^")
    catch e
        println("Warning $labelstr: $e")
    end

    return hcat(y, ym, ya)
end

function plotStatsPerClass(Results::FluidResultStruct, Trace::TraceStruct,
        QN::QueueNetworkStruct, Params::ParameterStruct; fignbr=1, n_cdf=1000)

    cdf_data = zeros(n_cdf, 10, sum(Params.C))

    # Extract mean using Littles Law
    trmLL = getAvgRT(Results.sol_smooth.u[end], Results.ps_smooth, QN, Params)

    figure(fignbr)
    clf()
    for k = 1:sum(Params.C)
        subplot(2, sum(Params.C), k)
        i0 = findfirst(Trace.q[k][:, 1] .> tspan[1])
        i1 = (x -> isnothing(x) ? size(Trace.q[k], 1) : x)(findfirst(Trace.q[k][:, 1] .>= tspan[2]))

        q_mavg, _ = movingAvg(Trace.q[k][i0:i1-1, 2], weights(diff(Trace.q[k][i0:i1, 1])), window=w)
        step(Trace.q[k][i0:i1, 1], Trace.q[k][i0:i1, 2], where="post", alpha=0.5, label="data")
        plot(Trace.q[k][i0+1:i1, 1], q_mavg, "C0", label="mAvg")
        plot(Results.sol_min.t, Results.qf_min[k], "C1", label="fluid")
        plot(Results.sol_smooth.t, Results.qf_smooth[k], "C2", label="fluid smooth")
        plot(Results.sol_opt.t, Results.qf_opt[k], "C3", label="fluid opt")
        title("queue length")
        xlim(tspan)
        if k == 1; legend(); end

        subplot(2, sum(Params.C), k + sum(Params.C))
        x = range(0, stop=quantile(Trace.tw[k], 0.99), length=n_cdf)
        i = findfirst(x-> k <= x, cumsum(Params.C))
        cdf_data[:, 1, k] = x
        cdf_data[:, 2:4, k] = tr_cdf_plot_data(x, Trace.tw[k], 0.95)
        cdf_data[:, 5:7, k] = tr_cdf_plot(x, Results.rt_cdf_min[k], 0.95, "C1", "fit CDF")
        cdf_data[:, 8:10, k] = tr_cdf_plot(x, Results.rt_cdf_smooth[k], 0.95, "C2", "fit CDF smooth")
       
        tr_cdf_plot(x, Results.rt_cdf_opt[k], 0.95, "C3", "fit CDF opt")
        cf4(t) = tr_cdf_qdist(t, (t, x) -> tr_cdf(t, Params.ph[k].T', x, Params.K[i], Params.ph[k].π), 
            Trace.qd[i])

        plot([trmLL[k], trmLL[k]], [0, 1], "k--", label="LL")
        if k == 1; legend(); end
    end

    return cdf_data
end

function plotCDFOverCN(Results::FluidResultStruct, Trace::TraceStruct,
        QN::QueueNetworkStruct, Params::ParameterStruct; fignbr=2, n_cdf=1000)

    # Extract mean RT using Littles Law from xs_smooth
    trmLLCN = getAvgRTCN(Results.sol_smooth.u[end], Results.ps_smooth, QN, Params)

    cdf_CN_data = zeros(n_cdf, 10)

    figure(fignbr)
    clf()

    x = range(0, stop=quantile(Trace.twCN, 0.975), length=n_cdf)

    cdf_CN_data[:, 1] = x
    cdf_CN_data[:, 2:4] = tr_cdf_plot_data(x, Trace.twCN, 0.95)
    cdf_CN_data[:, 5:7] = tr_cdf_plot(x, Results.rt_cdf_CN_min, 0.95, "C1", "tw min")
    cdf_CN_data[:, 8:10] = tr_cdf_plot(x, Results.rt_cdf_CN_smooth, 0.95, "C2", "tw smooth")
    tr_cdf_plot(x, Results.rt_cdf_CN_opt, 0.95, "C3", "tw opt")

    plot([trmLLCN, trmLLCN], [0, 1], "k--", label="LL")

    legend()

    return cdf_CN_data

end