module Sources

using SkyCoords

export SourceList, PointSource, AstroDiff, Atmospheric
abstract type Source end

"Point source"
struct PointSource <: Source
    "Source name"
    name::String
    "Source location in ICRS coordinates"
    coord::ICRSCoords
    "Redshift"
    z::Float64
end

abstract type DiffuseSource <: Source end

struct AstroDiff <: DiffuseSource
    name::String
end

struct Atmospheric <: DiffuseSource
    name::String
end

"Collection of sources to analyse"
@kwdef struct SourceList
    "List of point sources"
    PS::Vector{PointSource} = []
    "Astrophysical diffuse model"
    AstroDiff::Union{DiffuseSource, Bool} = false
    "Atmospheric model"
    Atmospheric::Union{DiffuseSource, Bool} = false
    #SourceList(PS, AstroDiff, Atmospheric) = PS isa Vector ? new(PS, AstroDiff, Atmospheric) : new([PS], AstroDiff, Atmospheric)
end

end