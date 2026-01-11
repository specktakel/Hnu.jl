module ROI

using SkyCoords

abstract type RegionOfInterest
end

struct CircularROI <: RegionOfInterest
    RA
    DEC
    center::ICRSCoords
    radius
end
CircularROI(RA, DEC, radius) = CircularROI(RA, DEC, ICRSCoords(RA, DEC), radius)
CircularROI(center::ICRSCoords, radius) = CircularROI(center.ra, center.dec, center, radius)

end

