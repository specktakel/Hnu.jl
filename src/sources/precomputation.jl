module ExposureIntegral

using Interpolations
using Integrals: IntegralProblem, solve, QuadGKJL

function calculateExposure(spectrum, aeff, params)
    function integrand(logx, p)
        params = merge(p, (E = 10^logx, log10x = logx, log10x0 = log10(p.x0)))
        return exp(aeff(logx, 0.)) * spectrum(params)
    end

    sol = solve(IntegralProblem(integrand, (2, 9), params), QuadGKJL());
    return sol
end

function calcExpGrid(spectrum, aeff, paramGrid, params)
    out = Vector{Float64}(undef, length(paramGrid))
    for i = 1:length(paramGrid)
        out[i] = calculateExposure(spectrum, aeff, merge((gamma=paramGrid[i],), params)).u
    end
    return out
end


end
