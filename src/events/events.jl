module Events

using DelimitedFiles
using SkyCoords
using Unitful

import ..ROI: CircularROI

const BASEPATH = "/Users/David/.icecube_data/20210126_PS-IC40-IC86_VII/icecube_10year_ps/events"

const IC40_PATH = joinpath(BASEPATH, "IC40_exp.csv")
const IC59_PATH = joinpath(BASEPATH, "IC59_exp.csv")
const IC79_PATH = joinpath(BASEPATH, "IC79_exp.csv")
const IC86_I_PATH = joinpath(BASEPATH, "IC86_I_exp.csv")
const IC86_II_PATH = joinpath(BASEPATH, "IC86_II_exp.csv")
const IC86_III_PATH = joinpath(BASEPATH, "IC86_III_exp.csv")
const IC86_IV_PATH = joinpath(BASEPATH, "IC86_IV_exp.csv")
const IC86_V_PATH = joinpath(BASEPATH, "IC86_V_exp.csv")
const IC86_VI_PATH = joinpath(BASEPATH, "IC86_VI_exp.csv")
const IC86_VII_PATH = joinpath(BASEPATH, "IC86_VII_exp.csv")

const IC40 = 1
const IC59 = 2
const IC79 = 3
const IC86_I = 4
const IC86_II = 5


mutable struct EventList
    mjd::Vector
    ra::Vector
    dec::Vector
    ang_err::Vector
    coords::Vector
    energy::Vector
    detector::Vector
    N::Int
end

EventList(mjd, ra, dec, ang_err, energy, detector) = EventList(mjd, ra, dec, ang_err, ICRSCoords.(ra, dec), energy, detector, length(energy))

function loadEvents(fname::String)
    if "40" occursin(fname)
        detector = IC40
    elseif "59" occursin(fname)
        detector = IC59
    elseif "79" occursin(fname)
        detector = IC79
    elseif "IC86_II" occursin(fname)
        detector = IC86_II
    elseif "IC86_I" occursin(fname)
        detector = IC86_I
    else
        println("unsupported detector")
    end
    events = readdlm(fname, Float64, comments=true, comment_char='#')
    detector_vector = Vector{Int}(detector, lenght(events[:, 2]))
    EventList(events[:, 1], events[:, 4]u"deg", events[:, 5]u"deg", events[:, 3]u"deg", events[:, 2], detector_vector)
end

function selectEvents!(events::EventList, roi::CircularROI)
end


function selectEvents!(events::EventList, mask)
    events.mjd = events.mjd[mask]
    events.ra = events.ra[mask]
    events.dec = events.dec[mask]
    events.ang_err = events.ang_err[mask]
    events.energy = events.energy[mask]
    events.coords = events.coords[mask]
    events.N = length(events.energy)
end

end