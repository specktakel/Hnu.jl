module EffectiveArea

using DelimitedFiles
using Interpolations

include("./detector.jl")
using .Detector

IC40_path = "IC40_effectiveArea.csv"
IC59_path = "IC59_effectiveArea.csv"
IC79_path = "IC79_effectiveArea.csv"
IC86_I_path = "IC86_I_effectiveArea.csv"
IC86_II_path = "IC86_II_effectiveArea.csv"


struct EffectiveArea
    season
    log10eBins    # log10(E/GeV)
    c_log10eBins
    sinDecBins   # dimension less
    c_sinDecBins
    area    # in m2
end



function load_aeff(season)
    if season == IC40
        path = joinpath(BasePath, IC40_path)
    elseif season == IC59
        path = joinpath(BasePath, IC59_path)
    elseif season == IC79
        path = joinpath(BasePath, IC79_path)
    elseif season == IC86_I
        path = joinpath(BasePath, IC86_I_path)
    elseif season == IC86_II
        path = joinpath(BasePath, IC86_II_path)
    else
        error("No other detector implemented")
    end

    data = readdlm(path, Float64, comments=true, comment_char='#')
    log10eBins = sort!(collect((Set(round.(data[:, 1:2], sigdigits=2)))))
    sinDecBins = sind.(sort!(collect((Set(round.(data[:, 3:4], sigdigits=2))))))
    N_energy = length(log10eBins) - 1
    N_dir = length(sinDecBins) - 1
    area = reshape(data[:, 5], (N_energy, N_dir)) * 1e-4   # unit:m2
    c_log10eBins = (log10eBins[1:end-1] + log10eBins[2:end]) / 2
    c_sinDecBins = (sinDecBins[1:end-1] + sinDecBins[2:end]) / 2
    c_log10eBins[1] = log10eBins[1]
    c_log10eBins[end] = log10eBins[end]
    c_sinDecBins[1] = sinDecBins[1]
    c_sinDecBins[end] = sinDecBins[end]
    EffectiveArea(season, log10eBins, c_log10eBins, sinDecBins, c_sinDecBins, area)
end

function construct_aeff_interpolation(aeff::EffectiveArea)
    area = copy(aeff.area)
    nonzero_min = minimum(area[area .> 0.])
    area[area.==0.] .= 1e-2 * nonzero_min
    linear_interpolation(
        (aeff.c_log10eBins, aeff.c_sinDecBins),
        area,
        extrapolation_bc = 1e-2 * nonzero_min
    )
end

function construct_aeff_log_interpolation(aeff::EffectiveArea)
    area = copy(aeff.area)
    nonzero_min = minimum(area[area .> 0.])
    area[area.==0.] .= 1e-2 * nonzero_min
    linear_interpolation(
        (aeff.c_log10eBins, aeff.c_sinDecBins),
        log.(area),
        extrapolation_bc = log(1e-2 * nonzero_min)
    )
end



end
