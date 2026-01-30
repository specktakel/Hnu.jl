module Detector

using DelimitedFiles
using Interpolations

struct EffectiveArea
    season
    log10eBins
    c_log10eBins
    sinDecBins
    c_sinDecBins
    areaGrid
    interp
    interp_log
end

const IC40 = 1
const IC59 = 2
const IC79 = 3
const IC86_I = 4
const IC86_II = 5

BasePath = joinpath(homedir(), ".icecube_data/20210126_PS-IC40-IC86_VII/icecube_10year_ps/irfs")
IC40_path = "IC40_effectiveArea.csv"
IC59_path = "IC59_effectiveArea.csv"
IC79_path = "IC79_effectiveArea.csv"
IC86_I_path = "IC86_I_effectiveArea.csv"
IC86_II_path = "IC86_II_effectiveArea.csv"


function loadEffectiveArea(season)
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
    area = reshape(data[:, 5], (N_energy, N_dir))
    c_log10eBins = (log10eBins[1:end-1] + log10eBins[2:end]) / 2
    c_sinDecBins = (sinDecBins[1:end-1] + sinDecBins[2:end]) / 2
    c_log10eBins[1] = log10eBins[1]
    c_log10eBins[end] = log10eBins[end]
    c_sinDecBins[1] = sinDecBins[1]
    c_sinDecBins[end] = sinDecBins[end]
    interp = linear_interpolation((c_log10eBins, c_sinDecBins), area)
    nonzero_min = minimum(area[area .> 0.])
    area[area.==0.] .= 1e-2 * nonzero_min
    interp_log = linear_interpolation((c_log10eBins, c_sinDecBins), log.(area))
    EffectiveArea(season, log10eBins, c_log10eBins, sinDecBins, c_sinDecBins, area, interp, interp_log)
end


struct EnergyResolution
    #season
    c_log10eRecoBins
    c_log10eTrueBins
    eres
    interp
end

function loadEnergyResolution()
    ereco = Vector(range(start=1.05, stop=9.00, step=0.01))
    path = joinpath(@__DIR__, "../../inputs/eres_ic86ii.dat")
    eres = readdlm(path)
    path = joinpath(@__DIR__, "../../inputs/eres_ic86ii_true_e_binc.dat")
    etrue = readdlm(path)[:, 1]
    etrue[1] = 2.0
    etrue[end] = 9.0
    interp = linear_interpolation((ereco, etrue), eres)
    EnergyResolution(ereco, etrue, eres, interp)
end

struct DetectorModel
    season
    eres
    aeff
end

end