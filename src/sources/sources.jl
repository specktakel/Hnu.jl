module Sources

using SkyCoords

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

end