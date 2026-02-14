module ExposureIntegral

using Interpolations
using Integrals: IntegralProblem, solve, QuadGKJL

function calculateExposure(spectrum, aeff, params, dec)
    function integrand(logx, p)
        params = merge(p, (E = 10^logx, log10E = logx, log10E0 = log10(p.E0)))
        return exp(aeff(logx, dec)) * spectrum(params)
    end

    sol = solve(IntegralProblem(integrand, (2, 9), params), QuadGKJL());
    return sol
end

function calcExpGrid(spectrum, aeff, paramGrid, params, dec)
    out = Vector{Float64}(undef, length(paramGrid))
    for i = 1:length(paramGrid)
        out[i] = calculateExposure(spectrum, aeff, merge((gamma=paramGrid[i],), params), dec).u
    end
    return out
end


end
