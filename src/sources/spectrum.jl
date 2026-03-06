module Spectrum

using Integrals: IntegralProblem, solve, QuadGKJL

function powerlaw(params::NamedTuple)
    N0 = params.norm
    x = params.E
    x0 = params.E0
    gamma = params.gamma
    return N0 * (x / x0)^(-gamma)
end

function powerlaw_logdomain(params::NamedTuple)
    N0 = params.norm
    log10x = params.log10E
    gamma = params.gamma
    log10x0 = params.log10E0
    x0 = params.E0
    return N0 * x0 * (10^(log10x - log10x0))^(-gamma + 1) * log(10.)
end

function powerlaw_energyflux(params::NamedTuple, Emin, Emax)
    function integrand(logx, p)
        params = merge(p, (log10E=logx,))
        return powerlaw_logdomain(params) * p.E0
    end
    p = merge(params, (log10E0=log10(params.E0), gamma=params.gamma + 1.))
    sol = solve(IntegralProblem(integrand, (log10(Emin), log10(Emax)), p), QuadGKJL())
    return sol.u
end

function powerlaw_numberflux(params::NamedTuple, Emin, Emax)
    function integrand(logx, p)
        params = merge(p, (log10E=logx,))
        return powerlaw_logdomain(params)
    end
    p = merge(params, (log10E0=log10(params.E0),))
    sol = solve(IntegralProblem(integrand, (log10(Emin), log10(Emax)), p), QuadGKJL())
    return sol.u
end

function powerlaw_calc_norm(Nex, T, gamma, exp_func)
    return Nex / T / exp_func(gamma)
end

end