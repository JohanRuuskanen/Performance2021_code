
function genSimSettings_static(;scaling_int=[0.8, 1.2])


    nf = 2
    Kf = [8, 8]
    mf = ones(nf*4) .* repeat(Kf', 4)[:] .* rand(Uniform(scaling_int...), nf*4)
    scvf = [5.0, 0.5, 5.0, 0.5, 5.0, 0.5, 5.0, 0.5] .* rand(Uniform(scaling_int...), nf*4)

    nb = 4
    Kb = [4, 4, 4, 4]
    mb = ones(nb*2) * 4 .* repeat(Kb', 2)[:] .* rand(Uniform(scaling_int...), nb*2)
    scvb = 10*ones(nb*2) .* rand(Uniform(scaling_int...), nb*2)

    N = 50

    μf = mf ./ repeat(Kf', 4)[:]
    μb = mb ./ repeat(Kb', 2)[:]

    ma = 10 # 10
    scva = 1.0

    md = 25 # 25
    scvd = 1.0

    ## Run simulation script for given data

    @assert ma > 0 && scva > 0
    @assert md > 0 && scvd > 0

    @assert length(Kf) == nf
    @assert length(mf) == length(scvf) == nf*4
    @assert all(mf .> 0) && all(scvf .> 0)

    @assert length(Kb) == nb
    @assert length(mb) == length(scvb) == nb*2
    @assert all(mb .> 0) && all(scvb .> 0)

    return SimSettingsStruct(N, ma, scva, md, scvd, nf, Kf, mf, scvf, nb, Kb, mb, scvb), [μf; μb]

end

function perturbArrivalsSimSettings(s::SimSettingsStruct; ma=0, md=0, scvd=0)
    ma_new = (ma == 0 ? s.ma : ma)
    md_new = (md == 0 ? s.md : md)
    scvd_new = (scvd == 0 ? s.scvd : scvd)

    return SimSettingsStruct(s.N, ma_new, s.scva, md_new, scvd_new, s.nf, 
        s.Kf, s.mf, s.scvf, s.nb, s.Kb, s.mb, s.scvb)
end

function perturbHScaleSimSettings(s::SimSettingsStruct, nf_ind::BitArray{1}, nb_ind::BitArray{1}, 
        p_smooth::Array{Float64, 1}; scaling_int=[0.8, 1.2])

    nf = sum(nf_ind)
    nb = sum(nb_ind)
    
    nf_ind_c = repeat(nf_ind', 4, 1)[:]
    nb_ind_c = repeat(nb_ind', 2, 1)[:]


    if nf == s.nf && nb == s.nb
        return s, p_smooth
    end

    if nf - s.nf == -1 && nb == s.nb
        Kf = s.Kf[nf_ind]
        mf = s.mf[nf_ind_c]
        scvf = s.scvf[nf_ind_c]

        return SimSettingsStruct(s.N, s.ma, s.scva, s.md, s.scvd, 
            nf, Kf, mf, scvf, s.nb, s.Kb, s.mb, s.scvb), p_smooth[[true; nf_ind; nb_ind]]
    elseif nf == s.nf && nb - s.nb == -1 
        Kb = s.Kb[nb_ind]
        mb = s.mb[nb_ind_c]
        scvb = s.scvb[nb_ind_c]

        return SimSettingsStruct(s.N, s.ma, s.scva, s.md, s.scvd, 
            s.nf, s.Kf, s.mf, s.scvf, nb, Kb, mb, scvb), p_smooth[[true; nf_ind; nb_ind]]
    elseif nf - s.nf == 1 && nb == s.nb
        Kf = [s.Kf; 8]
        mf = [s.mf; 8*ones(4).*rand(Uniform(scaling_int...), 4)]
        scvf = [s.scvf; [5.0, 0.5, 5.0, 0.5].*rand(Uniform(scaling_int...), 4)]

        p_smooth_new = zeros(1+nb+nf)
        p_smooth_new[1:1+s.nf] = p_smooth[1:1+s.nf]
        p_smooth_new[2+s.nf:1+nf] .= mean(p_smooth[2:1+s.nf])
        p_smooth_new[2+nf:end] = p_smooth[2+s.nf:end]

        return SimSettingsStruct(s.N, s.ma, s.scva, s.md, s.scvd, 
            nf, Kf, mf, scvf, s.nb, s.Kb, s.mb, s.scvb), p_smooth_new

    elseif nf == s.nf && nb - s.nb == 1
        Kb = [s.Kb; 4]
        mb = [s.mb; 16*ones(2).*rand(Uniform(scaling_int...), 2)]
        scvb = [s.scvb; 10*ones(2).*rand(Uniform(scaling_int...), 2)]

        p_smooth_new = zeros(1+nb+nf)
        p_smooth_new[1:1+s.nf+s.nb] = p_smooth[1:1+s.nf+s.nb]
        p_smooth_new[2+s.nf+s.nb:1+s.nf+nb] .= mean(p_smooth[2+s.nf:1+s.nf+s.nb])
        
        return SimSettingsStruct(s.N, s.ma, s.scva, s.md, s.scvd, 
            s.nf, s.Kf, s.mf, s.scvf, nb, Kb, mb, scvb), p_smooth_new

    else
        error("No such scaling has been implemented")
    end

end

function perturbVScaleSimSettings(s::SimSettingsStruct, target::Array{Int64, 1}, scale::Float64)
    @assert 1 <= target[1] <= 2

    if target[1] == 1
        @assert 1 <= target[2] <= s.nf

        a = ones(s.nf)
        a[target[2]] = scale
        a = repeat(a', 4, 1)[:]

        mf = s.mf .* a

        return SimSettingsStruct(s.N, s.ma, s.scva, s.md, s.scvd, s.nf, 
            s.Kf, mf, s.scvf, s.nb, s.Kb, s.mb, s.scvb)
    else 
        @assert 1 <= target[2] <= s.nb

        a = ones(s.nb)
        a[target[2]] = scale
        a = repeat(a', 2, 1)[:]

        mb = s.mb .* a

        return SimSettingsStruct(s.N, s.ma, s.scva, s.md, s.scvd, s.nf, 
            s.Kf, s.mf, s.scvf, s.nb, s.Kb, mb, s.scvb)
    end

end