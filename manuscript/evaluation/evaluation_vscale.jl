
basepath = "" # Enter path to the top folder
filepath = joinpath(basepath, "manuscript/evaluation")
datapath = joinpath(filepath, "results_vscale")

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

## Test 3: Vertical scaling

MCsims_vscale = 500

printstyled("\n===== RUNNING VERTICAL SCALING TEST =====\n",bold=true, color=:magenta)


err_vec_vscale = Array{ErrorStruct, 1}(undef, MCsims_vscale)
pes_vec_vscale = Array{SimSettingsStruct, 1}(undef, MCsims_vscale)

for k = 1:MCsims_vscale

    printstyled("## Running test $k / $MCsims_vscale ##\n",bold=true, color=:green)

    seed = seed_init*(k+1) + 1023414

    target = [rand([1,2]), 0]
    target[2] = (target[1] == 1 ? rand([1, 2]) : rand([1, 4]))
    scaling = rand() < 0.5 ? 1 / rand(Uniform(1, 5)) : rand(Uniform(1, 5))

    perturbedSimSettings = perturbVScaleSimSettings(simSettings, target, scaling)

    err_vec_vscale[k] = compareResults(runSim(perturbedSimSettings, CN_subset, 
        seed, p_smooth=p_smooth)...)
    pes_vec_vscale[k] = perturbedSimSettings

end

save(joinpath(filepath, "results_vscale.jld"), 
    "err_vec_vscale", err_vec_vscale, 
    "pes_vec_vscale", pes_vec_vscale)

## Test 3:   Plot results

bins = 0:0.05:1.0

figure(2)
clf()
subplot(1, 3, 1)
title("vertical scaling")
hist(getfield.(err_vec_vscale, :q_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :q_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :q_maxp_opt)[:], bins, alpha=0.5)
subplot(1, 3, 2)
hist(getfield.(err_vec_vscale, :tra_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_maxp_opt)[:], bins, alpha=0.5)
subplot(1, 3, 3)
hist(getfield.(err_vec_vscale, :tra_CN_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_CN_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_CN_opt)[:], bins, alpha=0.5)
