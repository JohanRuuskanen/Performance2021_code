
basepath = "" # Enter path to the top folder
filepath = joinpath(basepath, "manuscript/example2")
datapath = joinpath(basepath, "manuscript/example2/results/")

include(joinpath(basepath, "src/main.jl"))

if !(@isdefined matlab_session)
    matlab_session = MSession()
end


## Parameters

# Suffix to name files with
suffix = "1"

# Simulation seed
seed_init = 1

# Moving window
w = 1000

# Timespan of LINE sim
timespan = [0, 200000]

# Subset 
CN_subset = ["Class1" "Queue2";
             "Class1" "Queue3"]

## Run initial simulation to fit parameter

Random.seed!(seed_init)

foreach(rm, joinpath.(datapath, readdir(datapath)))

mxcall(matlab_session, :cd, 1, filepath)
mxcall(matlab_session, :simulateExample2, 0, suffix, datapath, 
    1/0.2, seed_init, timespan)


## Extract initial data and examine errors

QN, Params = extractParameters(suffix, CN_subset)
Trace = getTraceData(suffix, QN, Params)

tspan = (0.0, maximum(maximum.(Trace.ta)))

Results = getFluidResults(Trace.p_opt, QN, Params)
p_smooth = Trace.p_opt

Errors = compareResults(Results, Trace)

printErrors(Errors)

## Plot results

# Plot transients for selected sim
cdf_data = plotStatsPerClass(Results, Trace, QN, Params)

# Plot CDF over CN subset
cdf_CN_data = plotCDFOverCN(Results, Trace, QN, Params)

## Save CDF data

headers = ["X" hcat([["trCN$k" "trmCN$(k)" "traCN$(k)"] for k in ["data", "min", "smooth"]]...)]
outputfile = joinpath(datapath, "../example2_data_cdf.csv")
open(outputfile, "w") do f
    writedlm(f, headers, ",")
    writedlm(f, lowercase.(string.(cdf_CN_data[:, :, 1])), ",")
end

## Function for repeated evaluations

function runSim(μ1::Float64, seed; p_smooth=[])
    foreach(rm, joinpath.(datapath, readdir(datapath)))
    mxcall(matlab_session, :simulateExample2, 0, suffix, datapath, 1/μ1, seed, timespan)
    QN, Params = extractParameters(suffix, CN_subset)
    Trace = getTraceData(suffix, QN, Params, printout=true, warn=false)
    tspan = (0.0, maximum(maximum.(Trace.ta)))
    Results = getFluidResults(Trace.p_opt, QN, Params, p_smooth=p_smooth)
    return Results, Trace
end

## Scale rate of mu1

printstyled("\n===== RUNNING μ1 TEST =====\n",bold=true, color=:magenta)

mu_vec = 0.05:0.05:0.95
res_vec = Array{Any, 2}(undef, length(mu_vec), 2)

for (k, μ) in enumerate(mu_vec)
    printstyled("## Running test $k / $(length(mu_vec)) ##\n",bold=true, color=:green)
    res_vec[k, :] .= runSim(μ, seed_init*(k+1), p_smooth=p_smooth)
end

## Plot results

qm_data = hcat(getfield.(res_vec[:, 2], :qm)...)
qm_min = hcat(getfield.(res_vec[:, 1], :qm_min)...)
qm_smooth = hcat(getfield.(res_vec[:, 1], :qm_smooth)...)
qm_opt = hcat(getfield.(res_vec[:, 1], :qm_opt)...)

p_opt = hcat(getfield.(res_vec[:, 1], :p_opt)...)
p_smooth = hcat(getfield.(res_vec[:, 1], :p_smooth)...)

traCN_data = findQuantileData.(getfield.(res_vec[:, 2], :twCN), 0.95)
traCN_min = [fzero(t -> y(t) - 0.95, 1.0, order=0) 
                for y in getfield.(res_vec[:, 1], :rt_cdf_CN_min)] 
traCN_smooth = [fzero(t -> y(t) - 0.95, 1.0, order=0) 
                for y in getfield.(res_vec[:, 1], :rt_cdf_CN_smooth)] 
traCN_opt = [fzero(t -> y(t) - 0.95, 1.0, order=0) 
                for y in getfield.(res_vec[:, 1], :rt_cdf_CN_opt)] 


figure(3)
clf()

for k = 1:Params.Q
    subplot(Params.Q, 2, (k-1)*2 + 1)
    plot(mu_vec, qm_data[k, :], "C0", label="data")
    plot(mu_vec, qm_min[k, :], "C1", label="fluid")
    plot(mu_vec, qm_smooth[k, :], "C2", label="imp fluid")
    plot(mu_vec, qm_opt[k, :], "C3", label="imp fluid opt")
    title("Mean queue length")
    legend()

    if Params.Disc[k] == "PS"
        subplot(Params.Q, 2, (k-1)*2 + 2)
        plot(mu_vec, p_opt[k, :]', "C3", label="opt")
        plot(mu_vec, p_smooth[k, :]', "C2", label="p")
        title("p_opt")
        legend()
    end
end

figure(4)
plot(mu_vec, traCN_data, label="Data")
plot(mu_vec, traCN_min, label="fluid")
plot(mu_vec, traCN_smooth, label="imp fluid")
plot(mu_vec, traCN_opt, label="imp fluid opt")
title("95th percentile of response times over CN")

## Save data

headers = hcat(["X"], ["psb$k" for k = 1:Params.Q]..., 
    ["q$(k)data" for k = 1:Params.Q]...,
    ["q$(k)min" for k = 1:Params.Q]...,
    ["q$(k)smooth" for k = 1:Params.Q]...,
    ["q$(k)opt" for k = 1:Params.Q]...,
    "traCNdata", "traCNmin", "traCNsmooth", "traCNopt")

outputfile = joinpath(datapath, "../example2_data.csv")
open(outputfile, "w") do f
    writedlm(f, headers, ",")
    writedlm(f, vcat(mu_vec', p_opt, 
            qm_data, qm_min, qm_smooth, qm_opt,
            traCN_data', traCN_min', traCN_smooth', traCN_opt')', 
        ",")
end

