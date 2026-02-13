module LifeTime

using Unitful

include("./detector.jl")
using .Detector: IC40, IC59, IC79, IC86_I, IC86_II

const BasePath = joinpath(homedir(), ".icecube_data/20210126_PS-IC40-IC86_VII/icecube_10year_ps/uptime")

const IC40_path = "IC40_exp.csv"
const IC59_path = "IC59_exp.csv"
const IC79_path = "IC79_exp.csv"
const IC86_I_path = "IC86_I_exp.csv"
const IC86_II_path = "IC86_II_exp.csv"
const IC86_III_path = "IC86_III_exp.csv"
const IC86_IV_path = "IC86_IV_exp.csv"
const IC86_V_path = "IC86_V_exp.csv"
const IC86_VI_path = "IC86_VI_exp.csv"
const IC86_VII_path = "IC86_VII_exp.csv"

function loadUptime(season::Int)
    if season == IC40
        fname = joinpath(BasePath, IC40_path)
    elseif season == IC59
        fname = joinpath(BasePath, IC59_path)
    elseif season == IC79
        fname = joinpath(BasePath, IC79_path)
    elseif season == IC86_I
        fname = joinpath(BasePath, IC86_I_path)
    elseif season == IC86_II
        fname = joinpath(BasePath, IC86_II_path)
    else
        pritln("unsupported detector")
    end
    if season != IC86_II
        uptime = readdlm(fname, Float64, comments=true, comment_char='#')
    else
        uptime = readdlm(fname, Float64, comments=true, comment_char='#')
        for f in [IC86_III_path, IC86_IV_path, IC86_V_path, IC86_VI_path, IC86_VII_path]
            uptime = vcat(uptime, readdlm(joinpath(BasePath, f), Float64, comments=true, comment_char='#'))
        end
    end
end

function lifetime_from_dm(season)
    uptime = loadUptime(season)
    sum(uptime[:, 2] .- uptime[:, 1])u"d"
end

function mjd_from_dm(season)
    uptime = loadUptime(season)
    min = minimum(uptime[:, 1])
    max = maximum(uptime[:, 2])
    return (min, max)
end

function time_from_mjd(mjd_min, mjd_max)
    if mjd_max < mjd_min
        println("mjd_max < mjd_min")
        return 0.
    end
    if mjd_min >= maximum(uptime)
        println("mjd_min >= minimum(uptime)")
        return 0.
    elseif mjd_max <= minimum(uptime)
        println("mjd_max <= minimum(uptime)")
        return 0.
    end
    diffs = uptime[:, 2] .- uptime[:, 1]

    start_idx = searchsortedfirst(uptime[:, 1], mjd_min)
    end_idx = searchsortedfirst(uptime[:, 2], mjd_max)

    println(start_idx)
    println(end_idx)

    sum(diffs[start_idx:end_idx])
end

end
