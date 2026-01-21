module Detector

struct EffectiveArea
    season
    log10eBins
    c_log10eBins
    sinDecBins
    c_sinDecBins
    areaGrid
    interp
end

const IC40 = 1
const IC59 = 2
const IC79 = 3
const IC86_I = 4
const IC86_II = 5

BasePath = "/Users/David/.icecube_data/20210126_PS-IC40-IC86_VII/icecube_10year_ps/irfs"
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
    c_log10eBins = (ebins[1:end-1] + ebins[2:end]) / 2
    c_sinDecBins = (sinDecBins[1:end-1] + sinDecBins[2:end]) / 2
    c_log10eBins[1] = log10eBins[1]
    c_log10eBins[end] = log10eBins[end]
    c_sinDecBins[1] = sinDecBins[1]
    c_sinDecBins[end] = sinDecBins[end]
    interp = linear_interpolation((c_log10eBins, c_sinDecBins), area)
    EffectiveArea(season, log10eBins, c_log10eBins, sinDecBins, c_sinDecBins, area, interp)
end

end