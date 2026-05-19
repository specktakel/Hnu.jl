module EnergyResolution

using Interpolations
using DelimitedFiles

include("./detector.jl")
import .Detector

struct ERes{T<:Real}
    c_log10eRecoBins::Vector{T}
    c_log10eTrueBins::Vector{T}
    eres::Array{T, 2}
end

function load_eres()
    ereco = Vector(range(start=1.05, stop=9.00, step=0.01))
    path = joinpath(@__DIR__, "../../inputs/eres_ic86ii.dat")
    eres = readdlm(path)
    path = joinpath(@__DIR__, "../../inputs/eres_ic86ii_true_e_binc.dat")
    etrue = readdlm(path)[:, 1]
    etrue[1] = 2.0
    etrue[end] = 9.0
    ERes(ereco, etrue, eres)
end

function load_energy_llh(eres::ERes)
    ereco = eres.c_log10eRecoBins
    etrue = eres.c_log10eTrueBins
    grid = eres.eres
    interp = linear_interpolation((ereco, etrue), grid)
    return interp
end

end
