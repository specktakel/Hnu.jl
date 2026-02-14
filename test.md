---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.1
  kernelspec:
    display_name: Julia 1.12
    language: julia
    name: julia-1.12
---

```julia
using Pkg
Pkg.activate()
Pkg.develop(path="../Hnu")
```

```julia
using Revise
using Hnu
using DelimitedFiles
using Interpolations
using BAT
using SkyCoords
using Unitful
using Plots
using DensityInterface
using IntervalSets
using Distributions
using Integrals
using LaTeXStrings
using Roots
using Hnu.CosmoUnits
#using ForwardDiff
```

```julia
eres = Hnu.Detector.loadEnergyResolution();
energy_llh = Hnu.Detector.loadEnergyLikelihood(eres);
```

```julia
energy_llh(2.0, 4.0)
```

```julia
events = Hnu.Events.loadEvents(5)
```

```julia
ps = ICRSCoords(77.35u"deg", 5.7u"deg")
```

```julia
ps.dec
```

```julia
roi = Hnu.ROI.CircularROI(ps, 5u"deg")
```

```julia
Hnu.Events.selectEvents!(events, roi)
```

```julia
events.coords[1].dec
```

```julia
aeff = Hnu.Detector.loadEffectiveArea(5);
```

```julia
log_interp_aeff = Hnu.Detector.constructEffectiveAreaLogInterpolation(aeff);
```

```julia
abstract type Source end

struct PointSource <: Source
    name::String
    coord::ICRSCoords
    z::Float64
end

pointsource = PointSource("txs", ps, 0.3365)

ps_coll = (ps=pointsource, spectrum=Hnu.Spectrum.powerLawLogDomain)
```

```julia
gamma_grid = Vector(range(1.0, 4.0, step=0.05));
```

```julia
typeof(gamma_grid)
```

```julia
Revise.revise()
```

```julia
Hnu.ExposureIntegral.calculateExposure(ps_coll.spectrum, log_interp_aeff, (norm=1e-14, gamma=2.0, x0=1e5), ps.dec).u
```

```julia
exp_grid = Hnu.ExposureIntegral.calcExpGrid(ps_coll.spectrum, log_interp_aeff, gamma_grid, (norm = 1e-14, x0=1e5, ));
```

```julia
plt = plot(gamma_grid, exp_grid)
plot!(plt, yaxis=("exposure", :log10))
plot!(plt, xaxis=(L"$\gamma$"))
```

```julia
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
```

```julia
EnergyFluxPowerLaw((norm=1e-14, gamma=2., E0=1e5), 1e2, 1e9)
NumberFluxPowerLaw((norm=1e-14, gamma=2., E0=1e5), 1e2, 1e9)
```

```julia

```

```julia
function calculateExposure(spectrum, aeff, params, dec)
    function integrand(logx, p)
        params = merge(p, (E = 10^logx, log10x = logx, log10x0 = log10(p.x0)))
        return exp(aeff(logx, dec)) * spectrum(params)
    end

    sol = solve(IntegralProblem(integrand, (2, 9), params), QuadGKJL());
    return sol
end
```

```julia
typeof(ustrip(u"Mpc", Cosmology.luminosity_distance(2.)))
```

```julia
find_zero
```

```julia
uconvert(u"m", 1u"Mpc")
```

```julia
Om = 0.3
H0 = 70u"km/s"
```

```julia
function llh(Etrue, events, aeff, eres)
    llh = 0.
    logereco = events.energy
    coords = events.coords
    ra = [coords[i].ra for i in 1:events.N]
    dec = [coords[i].dec for i in 1:events.N]
    #Etrue = params.Etrue
    logEtrue = log10(Etrue)

    for i = 1:events.N
        llh += aeff.interp_log(logEtrue, sin(dec[i]))
        llh += eres.interp(logereco[i], logEtrue)
        llh += log(Etrue^(-2.19))
    end
    llh
end
```

```julia
Etrue = logrange(1e2, 1e9, length=1000)
likelihood = [llh(E, events, aeff, eres) for E in Etrue]
```

```julia
plt = plot(Etrue, likelihood)
plot!(xaxis=("log10E", :log10))
```

```julia
likelihood = let aeff = aeff, eres = eres, events = events
    logfuncdensity(function (params)
        function loglike(params, events, aeff, eres)
            llh = 0.
            logereco = events.energy
            coords = events.coords
            ra = [coords[i].ra for i in 1:events.N]
            dec = [coords[i].dec for i in 1:events.N]
            Etrue = params.Etrue
            logEtrue = log10(Etrue)

            for i = 1:events.N
                llh += aeff.interp_log(logEtrue, sin(dec[i]))
                llh += eres.interp(logereco[i], logEtrue)
                llh += log(Etrue^(-2.13))
            end
            llh
        end
        loglike(params, events, aeff, eres)
    end)
end
```

```julia
true_par_values = (Etrue=290e3,)
```

```julia
logdensityof(likelihood, true_par_values)
```

```julia
LogUniform
```

```julia
prior = distprod(Etrue=Uniform(1e2, 1e9))
```

```julia
posterior = PosteriorMeasure(likelihood, prior)
```

```julia
samples = bat_sample(posterior, TransformedMCMC(proposal = RandomWalk(), nsteps = 10^4, nchains = 2)).result
```

```julia
unshaped_samples, f_flatten = bat_transform(Vector, samples)
```

```julia
using LazyReports
lazyreport(samples)
```

```julia
histogram(log10.(samples.v.Etrue).-3)
```

```julia

```
