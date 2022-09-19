
basepath = "" # Enter path to the top folder
filepath = joinpath(basepath, "manuscript/example1")
datapath = joinpath(basepath, "manuscript/example1/results/")

include(joinpath(basepath, "src/main.jl"))

if !(@isdefined matlab_session)
    matlab_session = MSession()
end

if !isdir(datapath)
    mkdir(datapath)
end

## Parameters

# Suffix to name files with
suffix = "1"

# Simulation seed
seed_init = 1

# Moving window
w = 1000

# Timespan of LINE sim
timespan = [0, 500000] 

# Parameters for transient simulation
itrs = 10
tend = 10

# Subset 
CN_subset = ["Class1" "Queue"]

## Run initial simulation to fit parameter

Random.seed!(seed_init)

foreach(rm, joinpath.(datapath, readdir(datapath)))

mxcall(matlab_session, :cd, 1, filepath)

mxcall(matlab_session, :simulateExample1, 0, suffix, datapath, 
    1/0.6, seed_init, timespan)

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

headers = ["X" hcat([["tr1$k" "trm1$(k)" "tra1$(k)"] for k in ["data", "min", "smooth"]]...)]
outputfile = joinpath(datapath, "../example1_data_cdf.csv")
open(outputfile, "w") do f
    writedlm(f, headers, ",")
    writedlm(f, lowercase.(string.(cdf_data[:, :, 1])), ",")
end

## Function for repeated evaluations

function runSim(λ::Float64, seed; p_smooth=[])
    foreach(rm, joinpath.(datapath, readdir(datapath)))
    mxcall(matlab_session, :simulateExample1, 0, suffix, datapath, 1/λ, seed, timespan)
    QN, Params = extractParameters(suffix, CN_subset)
    Trace = getTraceData(suffix, QN, Params, printout=true, warn=false)
    tspan = (0.0, maximum(maximum.(Trace.ta)))
    Results = getFluidResults(Trace.p_opt, QN, Params, p_smooth=p_smooth)
    return Results, Trace, QN, Params
end

## Scale rate of arrivals

printstyled("\n===== RUNNING ARRIVAL TEST =====\n",bold=true, color=:magenta)

arrivals = 0.05:0.025:0.90
res_vec = Array{Any, 2}(undef, length(arrivals), 4)

for (k, λ) in enumerate(arrivals)
    printstyled("## Running test $k / $(length(arrivals)) ##\n",bold=true, color=:green)
    res_vec[k, :] .= runSim(λ, seed_init*(k+1), p_smooth=p_smooth)
end

## Run and plot experiments for transient values

arrivalIdx = [11, 27] # Using values 0.3, 0.7
itrs = 250
ts = [
    [0, 15],
    [0, 90]
]

for (i, k) in enumerate(arrivalIdx)
    mxcall(matlab_session, :simulateExample1_transVals, 0, datapath, "$i",
        1/arrivals[k], itrs, ts[i])
end

data = zeros(1001)
figure(3)
clf()
for (i, k) in enumerate(arrivalIdx)
    
    dt = (ts[i][2] - ts[i][1]) / 1000
    steps = ts[i][1]:dt:ts[i][2]
    qm_sim = getQM("transientVals_$i.mat", steps)
    fluid_min, fluid_smooth, _ = fluidSol(res_vec[k, 2], res_vec[k, 3], 
        res_vec[k, 4], ts[i], steps)
    
    subplot(length(arrivalIdx), 1, i)
    plot(steps, qm_sim[1])
    plot(steps, fluid_min[1])
    plot(steps, fluid_smooth[1])
    plot(ts[i], (res_vec[k, 2].qm ).* ones(2), "k--")

    data = hcat(data, collect(steps), qm_sim..., fluid_min..., fluid_smooth...)
end
data = data[:, 2:end]

## Save transient value data

headers = permutedims(vcat([["t_$i", "dataq1_$i", "fminq1_$i", "fsmoothq1_$i"] 
    for i = 1:length(arrivalIdx)]...))

outputfile = joinpath(datapath, "../example1_data_trans.csv")
open(outputfile, "w") do f
    writedlm(f, headers, ",")
    writedlm(f, data, ",")
end

## Plot results

qm_data = hcat(getfield.(res_vec[:, 2], :qm)...)
qm_min = hcat(getfield.(res_vec[:, 1], :qm_min)...)
qm_smooth = hcat(getfield.(res_vec[:, 1], :qm_smooth)...)
qm_opt = hcat(getfield.(res_vec[:, 1], :qm_opt)...)

p_opt = hcat(getfield.(res_vec[:, 1], :p_opt)...)

tra_data = hcat([findQuantileData.(x, 0.95) for x in getfield.(res_vec[:, 2], :tw)]...)
tra_min = hcat([ [fzero(t -> y(t) - 0.95, 1.0, order=0) for y in x] 
                for x in getfield.(res_vec[:, 1], :rt_cdf_min)]...)
tra_smooth = hcat([ [fzero(t -> y(t) - 0.95, 1.0, order=0) for y in x] 
                for x in getfield.(res_vec[:, 1], :rt_cdf_smooth)]...)
tra_opt = hcat([ [fzero(t -> y(t) - 0.95, 1.0, order=0) for y in x] 
                for x in getfield.(res_vec[:, 1], :rt_cdf_opt)]...)

figure(4)
clf()
subplot(3, 1, 1)
plot(arrivals, qm_data', label="data")
plot(arrivals, qm_min', label="fluid")
plot(arrivals, qm_smooth', label="imp fluid")
plot(arrivals, qm_opt', label="imp fluid opt")
title("Mean queue length")

subplot(3, 1, 2)
plot(arrivals, p_opt')
ylim([0.9, 1.1])
title("p_opt")

subplot(3, 1, 3)
plot(arrivals, tra_data', label="Data")
plot(arrivals, tra_min', label="fluid")
plot(arrivals, tra_smooth', label="imp fluid")
plot(arrivals, tra_opt', label="imp fluid opt")
title("95th percentile of response times")
xlabel("Arrival rate")

## Save data

headers = hcat(["X"], ["psb$k" for k = 1:Params.Q]..., 
    ["q$(k)data" for k = 1:Params.Q]...,
    ["q$(k)min" for k = 1:Params.Q]...,
    ["q$(k)smooth" for k = 1:Params.Q]...,
    ["q$(k)opt" for k = 1:Params.Q]...,
    ["tra$(k)data" for k = 1:Params.Q]...,
    ["tra$(k)min" for k = 1:Params.Q]...,
    ["tra$(k)smooth" for k = 1:Params.Q]...,
    ["tra$(k)opt" for k = 1:Params.Q]...)

outputfile = joinpath(datapath, "../example1_data.csv")
open(outputfile, "w") do f
    writedlm(f, headers, ",")
    writedlm(f, vcat(arrivals', p_opt, 
            qm_data, qm_min, qm_smooth, qm_opt,
            tra_data, tra_min, tra_smooth, tra_opt)', 
        ",")
end