module Utils
using SkyCoords

function calculate_separation(events, ref)
    sep = Vector{Float64}(undef, length(events.coords))
    for i in 1:length(events.coords)
        sep[i] = separation(events.coords[i], ref)
    end
    return sep
end

end