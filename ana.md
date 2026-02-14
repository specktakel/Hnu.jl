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
roi_mask = Vector{Bool}(undef, events.N);

for i = 1:events.N
    roi_mask[i] = separation(events.coords[i], roi.center) <= ustrip(u"rad", roi.radius)
end
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
length(mask)
```

```julia
bg = Hnu.Sources.loadBackgroundSource(Hnu.Detector.IC86_II)
```

```julia
bg_llh = bg.likelihood[roi_mask]
bg_llh = bg_llh[mask];
```

```julia
function signal_llh(params)

end
```
