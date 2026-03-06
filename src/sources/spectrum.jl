module Spectrum

using Integrals: IntegralProblem, solve, QuadGKJL

function powerLaw(params::NamedTuple)
    N0 = params.norm
    x = params.E
    x0 = params.E0
    gamma = params.gamma
    return N0 * (x / x0)^(-gamma)
end

function powerLawLogDomain(params::NamedTuple)
    N0 = params.norm
    log10x = params.log10E
    gamma = params.gamma
    log10x0 = params.log10E0
    x0 = params.E0
    return N0 * x0 * (10^(log10x - log10x0))^(-gamma + 1) * log(10.)
end

function EnergyFluxPowerLaw(params::NamedTuple, Emin, Emax)
    function integrand(logx, p)
        params = merge(p, (log10E=logx,))
        return powerLawLogDomain(params) * p.E0
    end
    p = merge(params, (log10E0=log10(params.E0), gamma=params.gamma + 1.))
    sol = solve(IntegralProblem(integrand, (log10(Emin), log10(Emax)), p), QuadGKJL())
    return sol.u
end

function NumberFluxPowerLaw(params::NamedTuple, Emin, Emax)
    function integrand(logx, p)
        params = merge(p, (log10E=logx,))
        return powerLawLogDomain(params)
    end
    p = merge(params, (log10E0=log10(params.E0),))
    sol = solve(IntegralProblem(integrand, (log10(Emin), log10(Emax)), p), QuadGKJL())
    return sol.u
end

function calcNorm(Nex, T, gamma, exp_func)
    return Nex / T / exp_func(gamma)
end

end