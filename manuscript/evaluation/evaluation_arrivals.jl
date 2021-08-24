
basepath = "" # Enter path to the top folder
filepath = joinpath(basepath, "manuscript/evaluation")
datapath = joinpath(filepath, "results_arrivals")

include(joinpath(basepath, "src/main.jl"))

if !(@isdefined matlab_session)
    matlab_session = MSession()
end
mxcall(matlab_session, :cd, 1, filepath)

## Parameters

# Suffix to name files with
suffix = "1"

# Simulation seed
seed_init = 1

# Moving window
w = 1000

# Timespan of LINE sim
timespan = [0, 500000]

# Subset 
CN_subset = ["ClosedClass3" "Frontend1";
         "ClosedClass4" "Frontend1";
         "ClosedClass3" "Backend1"; 
         "ClosedClass3" "Backend2";
         "ClosedClass3" "Backend3"; 
         "ClosedClass3" "Backend4"]

## Run initial simulation to fit parameter

Random.seed!(seed_init)
simSettings, μfd = genSimSettings_static()

foreach(rm, joinpath.(datapath, readdir(datapath)))

mxcall(matlab_session, :simulateTwoTierModel, 0, suffix, datapath, 
    simSettings, seed_init, timespan)

## Extract initial data and examine errors

QN, Params = extractParameters(suffix, CN_subset)
Trace = getTraceData(suffix, QN, Params)

tspan = (0.0, maximum(maximum.(Trace.ta)))

Results = getFluidResults(Trace.p_opt, QN, Params)
p_smooth = Trace.p_opt

Errors = compareResults(Results, Trace)
printErrors(Errors)

## Function for repeated evaluations

function runSim(s::SimSettingsStruct, subCN::Array{String,2}, seed; p_smooth=[])
    foreach(rm, joinpath.(datapath, readdir(datapath)))
    mxcall(matlab_session, :simulateTwoTierModel, 0, suffix, datapath, s, seed, timespan)
    QN, Params = extractParameters(suffix, subCN)
    Trace = getTraceData(suffix, QN, Params, printout=true, warn=false)
    tspan = (0.0, maximum(maximum.(Trace.ta)))
    Results = getFluidResults(Trace.p_opt, QN, Params, p_smooth=p_smooth)
    return Results, Trace
end

## Test 1: Scale rate of arrivals and client delays

MCsims_arrivals = 500

printstyled("\n===== RUNNING ARRIVAL TEST =====\n",bold=true, color=:magenta)

err_vec_arrivals = Array{ErrorStruct, 1}(undef, MCsims_arrivals)
pes_vec_arrivals = Array{SimSettingsStruct, 1}(undef, MCsims_arrivals)

for k in 1:MCsims_arrivals
    printstyled("## Running test $k / $MCsims_arrivals ##\n",bold=true, color=:green)
    seed = seed_init*(k+1)
    ma = rand(Uniform(2, 20))
    md = rand(Uniform(5, 50))
    scvd = rand(Uniform(0.1, 10.0))

    perturbedSimSettings = perturbArrivalsSimSettings(simSettings, 
        ma=ma, md=md, scvd=scvd)
    err_vec_arrivals[k] = compareResults(runSim(perturbedSimSettings, CN_subset, 
        seed, p_smooth=p_smooth)...)
    pes_vec_arrivals[k] = perturbedSimSettings
end

save(joinpath(filepath, "results_arrivals.jld"), 
    "err_vec_arrivals", err_vec_arrivals, 
    "pes_vec_arrivals", pes_vec_arrivals)

## Test 1:   Plot results

bins = 0:0.05:1.0

figure(1)
clf()
subplot(1, 3, 1)
title("workload change")
hist(getfield.(err_vec_arrivals, :q_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :q_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :q_meanp_opt)[:], bins, alpha=0.5)
subplot(1, 3, 2)
hist(getfield.(err_vec_arrivals, :tra_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_meanp_opt)[:], bins, alpha=0.5)
subplot(1, 3, 3)
hist(getfield.(err_vec_arrivals, :tra_CN_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_CN_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_CN_opt)[:], bins, alpha=0.5)
