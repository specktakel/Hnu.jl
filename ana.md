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
```

```julia
eres = Hnu.Detector.loadEnergyResolution();
energy_llh = Hnu.Detector.loadEnergyLikelihood(eres);
events = Hnu.Events.loadEvents(5)
ps = ICRSCoords(77.35u"deg", 5.7u"deg")
roi = Hnu.ROI.CircularROI(ps, 5u"deg")
aeff = Hnu.Detector.loadEffectiveArea(5);
log_interp_aeff = Hnu.Detector.constructEffectiveAreaLogInterpolation(aeff);

abstract type Source end

struct PointSource <: Source
    name::String
    coord::ICRSCoords
    z::Float64
end

pointsource = PointSource("txs", ps, 0.3365)


ps_coll = (ps=pointsource, spectrum=Hnu.Spectrum.powerLawLogDomain)

gamma_grid = Vector(range(1.0, 4.0, step=0.05));

function buildExpFunc(gamma_grid, exp_grid)
    interp = linear_interpolation(gamma_grid, log.(exp_grid))
    return x -> exp(interp(x))
end

exp_grid = Hnu.ExposureIntegral.calcExpGrid(ps_coll.spectrum, log_interp_aeff, gamma_grid, (norm = 1, E0=1e5, ), ps.dec);

exposure_function = buildExpFunc(gamma_grid, exp_grid)

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
Hnu.Events.selectEvents!(events, roi)
```

```julia
MJD_min=56917.
MJD_max=57113.


mask_min = events.mjd .>= MJD_min
mask_max = events.mjd .<= MJD_max

mask = mask_min .* mask_max;


```

```julia
Hnu.Events.selectEvents!(events, mask)
```

```julia
spatial_llh = Hnu.Events.calcSpatialLikelihood(events, pointsource.coord)
```

```julia
sinDecPS = sin(ustrip(u"rad", pointsource.coord.dec))
```

```julia
bg = Hnu.Sources.loadBackgroundSource(Hnu.Detector.IC86_II)
```

```julia
bg_llh = bg.likelihood[roi_mask]
bg_llh = bg_llh[mask];
```

```julia
all_events = Hnu.Events.loadEvents(5)
```

```julia
T = 15590030.
```

```julia
bg_norm = log(1 / log(1e9 / 1e2) / T / events.N)    # add this to bg_llh
```

```julia
function signal_llh(params)
    ### loglike = -Nex + sum_i (log(Aeff(E_i)) + log(spectrum(gamma, E_i)) + log(spatial(i)) 
    Nex = params.Nex
    gamma = params.gamma
    E = params.E
    llh = zeros(Float64, events.N)
    norm = calcNorm(Nex, T, gamma, exposure_function)
    for i = 1:events.N
        llh[i] += Hnu.Spectrum.powerLaw((norm = norm, E0=1e5, gamma=gamma, E=E[i]))
        llh[i] += spatial_llh[i]
        llh[i] += log_interp_aeff(log10(E[i]), sinDecPS)
        llh[i] += energy_llh(events.energy[i], log10(E[i]))
    end
    return llh
end

function background_llh(params)
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
background_llh((Nex_bg=250., E=events.energy .+1.))
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
signal_llh((gamma=2.0, Nex=10, E = 10 .^(events.energy .+ 1.)))
```

```julia

```
