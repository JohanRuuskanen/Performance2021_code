basepath = "" # Enter path to the top folder
filepath = joinpath(basepath, "manuscript/evaluation")

include(joinpath(basepath, "src/main.jl"))

## Examine load scale test

d = load(joinpath(filepath, "results_arrivals.jld"))

err_vec_arrivals = d["err_vec_arrivals"]
pes_vec_arrivals = d["err_vec_arrivals"]

bins = 0:0.01:1.0

figure(1)
clf()
subplot(2, 3, 1)
title("Workload test")
hist(getfield.(err_vec_arrivals, :q_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :q_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :q_meanp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 2)
hist(getfield.(err_vec_arrivals, :q_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :q_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :q_maxp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 4)
hist(getfield.(err_vec_arrivals, :tra_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_meanp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 5)
hist(getfield.(err_vec_arrivals, :tra_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_maxp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 6)
hist(getfield.(err_vec_arrivals, :tra_CN_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_CN_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_arrivals, :tra_CN_opt)[:], bins, alpha=0.5)


## Examine horizontal scale test

d = load(joinpath(filepath, "results_hscale.jld"))

err_vec_hscale_up = d["err_vec_hscale"]
pes_vec_hscale_up = d["pes_vec_hscale"]
err_vec_hscale_down = d["err_vec_hscale_down"]
pes_vec_hscale_down = d["pes_vec_hscale_down"]

printstyled("Downscale results\n",bold=true, color=:magenta)
for k = 1:length(err_vec_hscale_down)
    printstyled("## [nf, nb] = [$(pes_vec_hscale_down[k].nf), $(pes_vec_hscale_down[k].nb)] ##\n",
        bold=true, color=:green)
    printErrors(err_vec_hscale_down[k])
    println("")
end


err_vec_hscale = [err_vec_hscale_up; err_vec_hscale_down]
pes_vec_hscale = [pes_vec_hscale_up; pes_vec_hscale_down]

bins = 0:0.01:0.6

figure(2)
clf()
subplot(2, 3, 1)
title("Horzizontal upscale test")
hist(getfield.(err_vec_hscale, :q_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :q_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :q_meanp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 2)
hist(getfield.(err_vec_hscale, :q_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :q_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :q_maxp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 4)
hist(getfield.(err_vec_hscale, :tra_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_meanp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 5)
hist(getfield.(err_vec_hscale, :tra_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_maxp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 6)
hist(getfield.(err_vec_hscale, :tra_CN_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_CN_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_hscale, :tra_CN_opt)[:], bins, alpha=0.5)

headers = hcat(["X"], ["q_meanp_min"], ["q_meanp_smooth"], ["q_meanp_opt"],
    ["tra_meanp_min"], ["tra_meanp_smooth"], ["tra_meanp_opt"],
    ["tra_CN_min"], ["tra_CN_smooth"], ["tra_CN_opt"])

outputfile = joinpath(filepath, "evaluation_hscale.csv")
open(outputfile, "w") do f
    writedlm(f, headers, ",")
    writedlm(f, hcat(ones(506), 
        getfield.(err_vec_hscale, :q_meanp_min)[:],
        getfield.(err_vec_hscale, :q_meanp_smooth)[:],
        getfield.(err_vec_hscale, :q_meanp_opt)[:],
        getfield.(err_vec_hscale, :tra_meanp_min)[:],
        getfield.(err_vec_hscale, :tra_meanp_smooth)[:],
        getfield.(err_vec_hscale, :tra_meanp_opt)[:],
        getfield.(err_vec_hscale, :tra_CN_min)[:],
        getfield.(err_vec_hscale, :tra_CN_smooth)[:],
        getfield.(err_vec_hscale, :tra_CN_opt)[:]), ",")
end

## Examine vertical scale test

d = load(joinpath(filepath, "results_vscale.jld"))

err_vec_vscale = d["err_vec_vscale"]
pes_vec_vscale = d["pes_vec_vscale"]

bins = 0:0.01:0.7

figure(3)
clf()
subplot(2, 3, 1)
title("Vertical scaling test")
hist(getfield.(err_vec_vscale, :q_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :q_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :q_meanp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 2)
hist(getfield.(err_vec_vscale, :q_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :q_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :q_maxp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 4)
hist(getfield.(err_vec_vscale, :tra_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_meanp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 5)
hist(getfield.(err_vec_vscale, :tra_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_maxp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 6)
hist(getfield.(err_vec_vscale, :tra_CN_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_CN_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec_vscale, :tra_CN_opt)[:], bins, alpha=0.5)

headers = hcat(["X"], ["q_meanp_min"], ["q_meanp_smooth"], ["q_meanp_opt"],
    ["tra_meanp_min"], ["tra_meanp_smooth"], ["tra_meanp_opt"],
    ["tra_CN_min"], ["tra_CN_smooth"], ["tra_CN_opt"])

outputfile = joinpath(filepath, "evaluation_vscale.csv")
open(outputfile, "w") do f
    writedlm(f, headers, ",")
    writedlm(f, hcat(ones(500), 
        getfield.(err_vec_vscale, :q_meanp_min)[:],
        getfield.(err_vec_vscale, :q_meanp_smooth)[:],
        getfield.(err_vec_vscale, :q_meanp_opt)[:],
        getfield.(err_vec_vscale, :tra_meanp_min)[:],
        getfield.(err_vec_vscale, :tra_meanp_smooth)[:],
        getfield.(err_vec_vscale, :tra_meanp_opt)[:],
        getfield.(err_vec_vscale, :tra_CN_min)[:],
        getfield.(err_vec_vscale, :tra_CN_smooth)[:],
        getfield.(err_vec_vscale, :tra_CN_opt)[:]), ",")
end

## Examine combined tests

d = load(joinpath(filepath, "results_arrivals.jld"))

err_vec = d["err_vec_arrivals"]
pes_vec = d["err_vec_arrivals"]

d = load(joinpath(filepath, "results_hscale.jld"))

err_vec = [err_vec; d["err_vec_hscale"]]
pes_vec = [pes_vec; d["pes_vec_hscale"]]


d = load(joinpath(filepath, "results_vscale.jld"))

err_vec = [err_vec; d["err_vec_vscale"]]
pes_vec = [pes_vec; d["pes_vec_vscale"]]


bins = 0:0.01:1.0

figure(1)
clf()
subplot(2, 3, 1)
title("Combined tests")
hist(getfield.(err_vec, :q_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :q_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :q_meanp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 2)
hist(getfield.(err_vec, :q_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :q_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :q_maxp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 4)
hist(getfield.(err_vec, :tra_meanp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :tra_meanp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :tra_meanp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 5)
hist(getfield.(err_vec, :tra_maxp_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :tra_maxp_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :tra_maxp_opt)[:], bins, alpha=0.5)
subplot(2, 3, 6)
hist(getfield.(err_vec, :tra_CN_min)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :tra_CN_smooth)[:], bins, alpha=0.5)
hist(getfield.(err_vec, :tra_CN_opt)[:], bins, alpha=0.5)

headers = hcat(["X"], ["q_meanp_min"], ["q_meanp_smooth"], ["q_meanp_opt"],
    ["tra_meanp_min"], ["tra_meanp_smooth"], ["tra_meanp_opt"],
    ["q_maxp_min"], ["q_maxp_smooth"], ["q_maxp_opt"],
    ["tra_maxp_min"], ["tra_maxp_smooth"], ["tra_maxp_opt"],
    ["tra_CN_min"], ["tra_CN_smooth"], ["tra_CN_opt"])

outputfile = joinpath(filepath, "evaluation_all.csv")
open(outputfile, "w") do f
    writedlm(f, headers, ",")
    writedlm(f, hcat(ones(1500), 
        getfield.(err_vec, :q_meanp_min)[:],
        getfield.(err_vec, :q_meanp_smooth)[:],
        getfield.(err_vec, :q_meanp_opt)[:],
        getfield.(err_vec, :tra_meanp_min)[:],
        getfield.(err_vec, :tra_meanp_smooth)[:],
        getfield.(err_vec, :tra_meanp_opt)[:],
        getfield.(err_vec, :q_maxp_min)[:],
        getfield.(err_vec, :q_maxp_smooth)[:],
        getfield.(err_vec, :q_maxp_opt)[:],
        getfield.(err_vec, :tra_maxp_min)[:],
        getfield.(err_vec, :tra_maxp_smooth)[:],
        getfield.(err_vec, :tra_maxp_opt)[:],
        getfield.(err_vec, :tra_CN_min)[:],
        getfield.(err_vec, :tra_CN_smooth)[:],
        getfield.(err_vec, :tra_CN_opt)[:]), ",")
end