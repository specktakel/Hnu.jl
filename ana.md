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
Pkg.add("ProfileView")
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
using StatsFuns: logsumexp
using ArraysOfArrays
using Profile
```

```julia
eres = Hnu.EnergyResolution.load_eres()
```

```julia
eres = Hnu.EnergyResolution.load_eres();
energy_llh = Hnu.EnergyResolution.load_energy_llh(eres);
events = Hnu.Events.load_events(5)
ps = ICRSCoords(77.35u"deg", 5.7u"deg")
roi = Hnu.ROI.CircularROI(ps, 5u"deg")
aeff = Hnu.EffectiveArea.load_aeff(5);
log_interp_aeff = Hnu.EffectiveArea.construct_aeff_log_interpolation(aeff);

abstract type Source end

struct PointSource <: Source
    name::String
    coord::ICRSCoords
    z::Float64
end


struct IRF
    energy_llh
    aeff
    log_interp_aeff
end


pointsource = PointSource("txs", ps, 0.3365)


ps_coll = (ps=pointsource, spectrum=Hnu.Spectrum.powerlaw_logdomain)

gamma_grid = Vector(range(1.0, 4.0, step=0.05));


exp_grid = Hnu.ExposureIntegral.calculate_exposure_grid(ps_coll.spectrum, log_interp_aeff, gamma_grid, (norm = 1, E0=1e5, ), ps.dec);

exposure_function = Hnu.ExposureIntegral.build_1d_exposure_function(gamma_grid, exp_grid)

function calcNorm(Nex, T, gamma, exp_func)
    return Nex / T / exp_func(gamma)
end
```

```julia
ustrip.(u"rad", events.ang_err).^2
```

```julia
roi_mask = Vector{Bool}(undef, events.N);

for i = 1:events.N
    roi_mask[i] = separation(events.coords[i], roi.center) <= ustrip(u"rad", roi.radius)
end
```

```julia
events.N
```

```julia
Hnu.Events.select_events!(events, roi)
```

```julia
MJD_min = 58010
MJD_max = 58020


mask_min = events.mjd .>= MJD_min
mask_max = events.mjd .<= MJD_max

mask = mask_min .* mask_max;


```

```julia
Hnu.Events.select_events!(events, mask)
```

```julia
spatial_llh = Hnu.Events.calc_spatial_llh(events, pointsource.coord)
```

```julia
sinDecPS = sin(ustrip(u"rad", pointsource.coord.dec))
```

```julia
bg = Hnu.Sources.load_background_source(Hnu.Detector.IC86_II)
```

```julia
bg_llh = bg.likelihood[roi_mask]
bg_llh = bg_llh[mask];
```

```julia
all_events = Hnu.Events.load_events(5)
```

```julia
# T = 15590030.
T = 848956.89
```

```julia
events.N
```

```julia
bg_norm = log(1 / log(1e9 / 1e2)) .+ log(761162) .- log(1.8998668e8) .- log(events.N)
```

```julia
"""
function signal_llh(params)
    ### loglike = -Nex + sum_i (log(Aeff(E_i)) + log(spectrum(gamma, E_i)) + log(spatial(i)) 
    Nex = params.Nex
    gamma = params.gamma
    E = params.E
    llh = zeros(Float64, events.N)
    norm = calcNorm(Nex, T, gamma, exposure_function)
    for i = 1:events.N
        llh[i] += log(Hnu.Spectrum.powerLaw((norm = norm, E0=1e5, gamma=gamma, E=E[i])))
        llh[i] += spatial_llh[i]
        llh[i] += log_interp_aeff(log10(E[i]), sinDecPS)
        llh[i] += energy_llh(events.energy[i], log10(E[i]))
    end
    return llh
end
"""

function build_signal_llh(T, exposure_function, spectrum, spatial_llh, log_interp_aeff, energy_llh)
    function signal_llh(params::NamedTuple)
        Nex = params.Nex
        gamma = params.gamma
        E = params.E
        llh = zeros(Float64, events.N)
        norm = calcNorm(Nex, T, gamma, exposure_function)
        for i = 1:events.N
            llh[i] += log(Hnu.Spectrum.powerlaw((norm = norm, E0=1e5, gamma=gamma, E=E[i])))
            llh[i] += spatial_llh[i]
            llh[i] += log_interp_aeff(log10(E[i]), sinDecPS)
            llh[i] += energy_llh(events.energy[i], log10(E[i]))
        end
        return llh
    end
    return x -> signal_llh(x)
end

function background_llh(params::NamedTuple)
    Nex_bg = params.Nex_bg
    E = params.E
    log_Nex_bg = log(Nex_bg)
    llh = zeros(Float64, events.N)
    for i = 1:events.N
        llh[i] = bg_llh[i] + bg_norm - log(E[i]) + log_Nex_bg
    end
    return llh
end
```

```julia
signal_llh = build_signal_llh(T, exposure_function, ps_coll.spectrum, spatial_llh, log_interp_aeff, energy_llh)
```

```julia
signal_llh((Nex=1., E=fill(1e3, events.N), Nex_bg=10., gamma=2.0))
```

```julia
function llh(params)
    event_llh = zeros(Float64, events.N, 2)
    event_llh[:, 1] = signal_llh(params)
    event_llh[:, 2] = background_llh(params)
    return event_llh
end
```

```julia
params = (Nex=1., Nex_bg=10., E=fill(1e3, events.N), gamma=2.0)
```

```julia
@profview llh(params)
```

```julia
Profile.print(format=:flat)
```

```julia
likelihood = let signal_llh = signal_llh, background_llh = background_llh
    logfuncdensity(function (params)
        function loglike(params)
            signal = signal_llh(params)
            bg = background_llh(params)

            llh = -params.Nex - params.Nex_bg
            for i=1:length(signal)
                llh += logsumexp(signal[i], bg[i])
            end
            llh
        end
        loglike(params)
    end)
end
```

```julia
prior = distprod(
    Nex = Uniform(0, 100),
    Nex_bg = Uniform(0, 100),
    gamma = truncated(Normal(2.3, 1.5), 1., 4.),
    E = fill(Uniform(1e2, 1e9), events.N),
)
```

```julia
posterior = PosteriorMeasure(likelihood, prior)
```

```julia
sampled = bat_sample(posterior, TransformedMCMC(proposal = RandomWalk(), nsteps = 10^4, nchains = 2))
```

```julia
sampled.evaluated.measure.likelihood
```

```julia
using LazyReports
lazyreport(samples)
```

```julia
? bat_sample
```

```julia
plot(
    samples, :(gamma), mean=true, std=true,
)
```

```julia
copied = samples.v.E[:, 1][:, 1, 1]
```

```julia
copied
```

```julia
energies = flatview(samples.v.E[:, 1])
```

```julia
plot(log10.(energies[10, :]), mean=true, nbins=50)
```

```julia
function test(a::Unitful.Energy)
    return a * 2
end

test(23u"GeV")
```

```julia
log10(ustrip(u"GeV", 123u"MeV"))
```

```julia
function outer(a)
    function inner(x)
        return x + a
    end
    return x -> inner(x)
end
```

```julia
func = outer(2)
```

```julia
func(2.2)
```

```julia

```
