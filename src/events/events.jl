module Events

using DelimitedFiles
using SkyCoords
using Unitful

using ..utils.ROI

const BASEPATH = "/Users/David/.icecube_data/20210126_PS-IC40-IC86_VII/icecube_10year_ps/events"
const IC40 = string(BASEPATH, "/", "IC40_exp.csv")

mutable struct EventList
    mjd::Vector
    ra::Vector
    dec::Vector
    ang_err::Vector
    coords::Vector
    energy::Vector
    N::Int
end

EventList(mjd, ra, dec, ang_err, energy) = EventList(mjd, ra, dec, ang_err, ICRSCoords.(ra, dec), energy, length(energy))

function loadEvents(fname::String)
    events = readdlm(fname, Float64, comments=true, comment_char='#')
    EventList(events[:, 1], events[:, 4]u"deg", events[:, 5]u"deg", events[:, 3]u"deg", events[:, 2])
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