module Sources

using SkyCoords
using DelimitedFiles

#include("./spectrum.jl")
#import .Spectrum
include("../detector/detector.jl")
import .Detector: IC86_II

export SourceList, PointSource, AstroDiff, Atmospheric

abstract type Source end

struct PointSource <: Source
    name::String
    coord::ICRSCoords
    z::Float64
end

abstract type DiffuseSource <: Source end

struct AstroDiff <: DiffuseSource
    name::String
end

struct Atmospheric <: DiffuseSource
    name::String
end

@kwdef struct SourceList
    PS::Vector{PointSource} = []
    AstroDiff::Union{DiffuseSource, Bool} = false
    Atmospheric::Union{DiffuseSource, Bool} = false
    #SourceList(PS, AstroDiff, Atmospheric) = PS isa Vector ? new(PS, AstroDiff, Atmospheric) : new([PS], AstroDiff, Atmospheric)
end

struct BackgroundSource <: Source
    season
    likelihood
end

function load_background_source(season::Int)
    if season != IC86_II
        println("season not implemented")
        return 0.
    end

    llh = readdlm(joinpath(@__DIR__, "../../inputs/ic86_ii_bg_llh.csv"), Float64, comments=false)
    bg = BackgroundSource(season, llh)
    return bg
end

end