
basepath = "" # Enter path to the top folder
filepath = joinpath(basepath, "manuscript/evaluation")
datapath = joinpath(filepath, "results_hscale")

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

## Test 2: Horizontal scaling

function getNewSubset(nb)
    fs = ["ClosedClass3" "Frontend1";
        "ClosedClass4" "Frontend1"]
    return vcat(fs, [["ClosedClass3" "Backend$k"] for k = 1:nb]...)
end

MCsims_hscale = 500

printstyled("\n===== RUNNING DOWNSCALING TEST =====\n",bold=true, color=:magenta)

err_vec_hscale_down = Array{ErrorStruct, 1}(undef, 6)
pes_vec_hscale_down = Array{SimSettingsStruct, 1}(undef, 6)

for k = 1:6
    printstyled("## Running Downscale $k / 6 ##\n",bold=true, color=:green)

    seed = seed_init*(k+1) + 12345

    nf_ind = Bool.([1, 1])
    nb_ind = Bool.([1, 1, 1, 1])
    if k <= 2
        nf_ind[k] = 0
    else
        nb_ind[k-2] = 0
    end

    perturbedSimSettings, p_smooth_pert = perturbHScaleSimSettings(simSettings, 
        nf_ind, nb_ind, p_smooth)

    err_vec_hscale_down[k] = compareResults(runSim(perturbedSimSettings, 
        getNewSubset(sum(nb_ind)), seed, p_smooth=p_smooth_pert)...)
    pes_vec_hscale_down[k] = perturbedSimSettings
end

printstyled("\n===== RUNNING UPSCALING TEST =====\n",bold=true, color=:magenta)

err_vec_hscale = Array{ErrorStruct, 1}(undef, MCsims_hscale)
pes_vec_hscale = Array{SimSettingsStruct, 1}(undef, MCsims_hscale)

for k = 1:MCsims_hscale
    printstyled("## Running Upscale $k / $MCsims_hscale ##\n",bold=true, color=:green)

    seed = seed_init*(k+1) + 54323

    nf_ind = Bool.([1, 1])
    nb_ind = Bool.([1, 1, 1, 1])
    if rand([1, 2]) == 1
        nf_ind = Bool.([1, 1, 1])
    else
        nb_ind = Bool.([1, 1, 1, 1, 1])
    end

    perturbedSimSettings, p_smooth_pert = perturbHScaleSimSettings(simSettings, 
        nf_ind, nb_ind, p_smooth, scaling_int=[0.5, 1.5])

    err_vec_hscale[k] = compareResults(runSim(perturbedSimSettings, 
        getNewSubset(sum(nb_ind)), seed, p_smooth=p_smooth_pert)...)
    pes_vec_hscale[k] = perturbedSimSettings

end

save(joinpath(filepath, "results_hscale.jld"), 
    "err_vec_hscale", err_vec_hscale, 
    "pes_vec_hscale", pes_vec_hscale,
    "err_vec_hscale_down", err_vec_hscale_down, 
    "pes_vec_hscale_down", pes_vec_hscale_down)

## Test 2:   Print and plot results

printstyled("Downscale results\n",bold=true, color=:magenta)
for k = 1:length(err_vec_hscale_down)
    printstyled("## [nf, nb] = [$(pes_vec_hscale_down[k].nf), $(pes_vec_hscale_down[k].nb)] ##\n",
        bold=true, color=:green)
    printErrors(err_vec_hscale_down[k])
    println("")
end

bins = 0:0.05:1.0

figure(3)
clf()
subplot(1, 3, 1)
title("horisontal upscaling")
hist(getfield.(err_vec_hscale, :q_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :q_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :q_maxp_opt)[:], bins, alpha=0.5)
subplot(1, 3, 2)
hist(getfield.(err_vec_hscale, :tra_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_maxp_opt)[:], bins, alpha=0.5)
subplot(1, 3, 3)
hist(getfield.(err_vec_hscale, :tra_CN_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_CN_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_CN_opt)[:], bins, alpha=0.5)
