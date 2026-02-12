module Spectrum

function powerLaw(params::NamedTuple)
    N0 = params.norm
    x = params.E
    x0 = params.E0
    gamma = params.gamma
    return N0 * (x / x0)^(-gamma)
end

function powerLawLogDomain(params::NamedTuple)
    N0 = params.norm
    log10x = params.log10x
    gamma = params.gamma
    log10x0 = params.log10x0
    x0 = params.x0
    return N0 * x0 * (10^(log10x - log10x0))^(-gamma + 1) * log(10.)
end

end