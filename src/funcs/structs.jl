
struct ParameterStruct
    Queues::Array{String, 1}
    Classes::Array{String, 1}
    Disc::Array{String, 1}
    Q::Int64
    C::Array{Int64, 1}
    C_idx::Array{UnitRange{Int64},1}
    C_set::Array{Array{String, 1}, 1}
    K::Array{Float64,1}
    subCN::Array{String, 2}
    subCNa::Array{Int64, 1}
    subCNaC::Array{Int64, 1}
    S::Array{Int64, 1}
    M::Array{Int64, 1}
    Mc::Array{Int64, 1}
    N::Array{Int64, 1}
    ph::Array{EMpht.PhaseType, 1}
    ph_arrivals::Array{EMpht.PhaseType, 1}
end

struct QueueNetworkStruct
    Φ::AbstractArray{Float64, 2}
    A::AbstractArray{Float64, 2}
    B::AbstractArray{Float64, 2}
    P::AbstractArray{Float64, 2}
    W::AbstractArray{Float64, 2}
    E::AbstractArray{Float64, 2}
    λ::AbstractArray{Float64, 2}
    Pb::AbstractArray{Float64, 2}
    Wb::AbstractArray{Float64, 2}
    X0::Array{Float64, 1}
end

struct TraceStruct
    ta::Array{Array{Float64, 1}}
    td::Array{Array{Float64, 1}}
    taCN::Array{Float64, 1}
    tdCN::Array{Float64, 1}
    q::Array{Array{Float64, 2}}
    q_node::Array{Array{Float64, 2}}
    qd::Array{Array{Float64, 1}}
    qCN::Array{Float64, 2}
    tw::Array{Array{Float64, 1}}
    twCN::Array{Float64, 1}
    qm::Array{Float64, 1}
    qm_node::Array{Float64, 1}
    p_opt::Array{Float64, 1}
    ps::Array{Float64, 1}
end

struct FluidResultStruct
    sol_min::ODESolution
    qf_min::Array{Array{Float64, 1}}
    qm_min::Array{Float64, 1}
    ps_min::Array{Float64, 1}
    rt_cdf_min::Array{Function, 1}
    rt_cdf_CN_min::Function

    sol_smooth::ODESolution
    qf_smooth::Array{Array{Float64, 1}}
    qm_smooth::Array{Float64, 1}
    ps_smooth::Array{Float64, 1}
    p_smooth::Array{Float64, 1}
    rt_cdf_smooth::Array{Function, 1}
    rt_cdf_CN_smooth::Function

    sol_opt::ODESolution
    qf_opt::Array{Array{Float64, 1}}
    qm_opt::Array{Float64, 1}
    ps_opt::Array{Float64, 1}
    p_opt::Array{Float64, 1}
    rt_cdf_opt::Array{Function, 1}
    rt_cdf_CN_opt::Function
end

struct ErrorStruct
    q_meanp_min::Float64
    q_meanp_smooth::Float64
    q_meanp_opt::Float64

    q_maxp_min::Float64
    q_maxp_smooth::Float64
    q_maxp_opt::Float64

    trm_meanp_min::Float64
    trm_meanp_smooth::Float64
    trm_meanp_opt::Float64

    trm_maxp_min::Float64
    trm_maxp_smooth::Float64
    trm_maxp_opt::Float64

    tra_meanp_min::Float64
    tra_meanp_smooth::Float64
    tra_meanp_opt::Float64

    tra_maxp_min::Float64
    tra_maxp_smooth::Float64
    tra_maxp_opt::Float64

    trm_CN_min::Float64
    trm_CN_smooth::Float64
    trm_CN_opt::Float64

    tra_CN_min::Float64
    tra_CN_smooth::Float64
    tra_CN_opt::Float64
end

struct SimSettingsStruct
    N::Int64

    ma::Float64
    scva::Float64

    md::Float64
    scvd::Float64
    
    nf::Int64
    Kf::Array{Int64, 1}
    mf::Array{Float64, 1}
    scvf::Array{Float64, 1}

    nb::Int64
    Kb::Array{Int64, 1}
    mb::Array{Float64, 1}
    scvb::Array{Float64, 1}
end