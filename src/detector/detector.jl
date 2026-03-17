module Detector

using DelimitedFiles
using Interpolations

export IC40, IC59, IC79, IC86_I, IC86_II, DetectorModel, BasePath

const IC40 = 1
const IC59 = 2
const IC79 = 3
const IC86_I = 4
const IC86_II = 5

BasePath = joinpath(homedir(), ".icecube_data/20210126_PS-IC40-IC86_VII/icecube_10year_ps/irfs")


struct DetectorModel
    season
    eres
    aeff
end

end