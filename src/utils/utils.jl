module Utils
using SkyCoords

function calculate_separation(events, ref)
    sep = zeros(length(events.coords))
    idxs = eachindex(events.coords)
    for i in idxs
        sep[i] = separation(events.coords[i], ref)
    end
    return sep
end

end