module Events

using DelimitedFiles
using SkyCoords
using Unitful

import ..ROI: CircularROI

const BASEPATH = joinpath(homedir(), ".icecube_data/20210126_PS-IC40-IC86_VII/icecube_10year_ps/events")

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

"""
    loadEvent(season)

Load all events recorded during the provided season.
"""
function loadEvents(season)
    fname = IC40_PATH
    if season == IC40
        fname = IC40_PATH
    elseif season == IC59
        fname = IC59_PATH
    elseif season == IC79
        fname = IC79_PATH
    elseif season == IC86_I
        fname = IC86_I_PATH
    elseif season == IC86_II
        fname == IC86_II_PATH
    else
        println("unsupported detector")
    end
    if season != IC86_II
        events = readdlm(fname, Float64, comments=true, comment_char='#')
    else
        events = readdlm(IC86_II_PATH, Float64, comments=true, comment_char='#') 
        for fname in [IC86_III_PATH, IC86_IV_PATH, IC86_V_PATH, IC86_VI_PATH, IC86_VII_PATH]
            events = vcat(events, readdlm(fname, Float64, comments=true, comment_char='#'))
        end
    end
    detector_vector = fill(season, length(events[:, 2]))
    EventList(events[:, 1], events[:, 4]u"deg", events[:, 5]u"deg", events[:, 3]u"deg", events[:, 2], detector_vector)
end

"""
    selectEvents!(events::EventList, roi::CircularROI)

Restrict events to the provided roi.
"""
function selectEvents!(events::EventList, roi::CircularROI)
    mask = Vector{Bool}(undef, events.N);

    for i = 1:events.N
        mask[i] = separation(events.coords[i], roi.center) <= ustrip(u"rad", roi.radius)
    end
    selectEvents!(events, mask)
end


"""
    selectEvents!(events::EventList, mask)

Select events based on provided mask
"""
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