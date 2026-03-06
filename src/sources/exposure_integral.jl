module ExposureIntegral

using Interpolations
using Integrals: IntegralProblem, solve, QuadGKJL
using Unitful

function calculateExposure(spectrum, aeff, params, dec, Emin::Unitful.Energy=1e2u"GeV", Emax::Unitful.Energy=1e9u"GeV")
    function integrand(logx, p)
        params = merge(p, (E = 10^logx, log10E = logx, log10E0 = log10(p.E0)))
        return exp(aeff(logx, dec)) * spectrum(params)
    end

    sol = solve(IntegralProblem(integrand, (log10(ustrip(u"GeV", Emin)), log10(ustrip(u"GeV", Emax))), params), QuadGKJL());
    return sol
end

function calcExpGrid(spectrum, aeff, paramGrid, params, dec)
    out = Vector{Float64}(undef, length(paramGrid))
    for i = 1:length(paramGrid)
        out[i] = calculateExposure(spectrum, aeff, merge((gamma=paramGrid[i],), params), dec).u
    end
    return out
end

function build1DExpFunc(param_grid, exp_grid)
    interp = linear_interpolation(param_grid, log.(exp_grid))
    return x -> exp(interp(x))
end

end
