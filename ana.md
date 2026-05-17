---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: Julia 1.12
    language: julia
    name: julia-1.12
---

```julia
using Pkg
Pkg.activate(".")
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
using AdvancedHMC
using ADTypes, ForwardDiff
```

```julia
# 10 events simulated using hierarchical_nu, index is 2.1, sky position of TXS
N = 36
data = readdlm("ps_2.1_10_bg_26_events.csv", Float64, comments=true, comment_char='#');
```

```julia
events = Hnu.Events.EventList(fill(99, N), data[:, 3]u"deg", data[:, 4]u"deg", data[:, 2]u"deg", data[:, 1], fill(5, N));
```

```julia
eres = Hnu.EnergyResolution.load_eres();
energy_llh = Hnu.EnergyResolution.load_energy_llh(eres);
# events = Hnu.Events.load_events(5)
ps = ICRSCoords(77.35u"deg", 5.7u"deg")
roi = Hnu.ROI.CircularROI(ps, 5u"deg")
aeff = Hnu.EffectiveArea.load_aeff(5);
log_interp_aeff = Hnu.EffectiveArea.construct_aeff_log_interpolation(aeff);

pointsource = Hnu.Sources.PointSource("txs", ps, 0.3365)

ps_coll = (ps=pointsource, spectrum=Hnu.Spectrum.powerlaw)

gamma = Hnu.Parameters.LinearParameter("index", 1.0, 4.0, 0.05)
gamma_grid = Hnu.Parameters.make_parameter_grid(gamma)

exp_grid = Hnu.ExposureIntegral.calculate_exposure_grid(Hnu.Spectrum.powerlaw_logdomain, log_interp_aeff, gamma_grid, (norm = 1, E0=1e5, ), ps.dec);

exposure_function = Hnu.ExposureIntegral.build_1d_exposure_function(gamma_grid, exp_grid)

function make_calcNorm_function(T, exp_func)
    function calcNorm(Nex, gamma)
        return Nex / T / exp_func(gamma)
    end

    return calcNorm
end

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
"""
MJD_min = 58010
MJD_max = 58020


mask_min = events.mjd .>= MJD_min
mask_max = events.mjd .<= MJD_max

mask = mask_min .* mask_max;
Hnu.Events.select_events!(events, mask)
"""
```

```julia
spatial_llh = Hnu.Events.calc_spatial_llh(events, pointsource.coord);
```

```julia
sinDecPS = sin(ustrip(u"rad", pointsource.coord.dec))
```

```julia
#bg = Hnu.Sources.load_background_source(Hnu.Detector.IC86_II)
#bg_llh = bg.likelihood[roi_mask]
#bg_llh = bg_llh[mask];
```

```julia
bg_llh = log.(readdlm("ps_2.1_10_bg_26_events_bg_likelihood.csv", comments=true, comment_char='#')[:, 1]);
```

```julia
#all_events = Hnu.Events.load_events(5)
```

```julia
# T = 15590030.
# T = 848956.89
T = 1829442.4

calcNorm = make_calcNorm_function(T, exposure_function)
```

```julia
bg_norm = log(1 / log(1e9 / 1e2)) .+ log(761162) .- log(1.8998668e8) .- log(events.N)
```

```julia
function build_signal_llh(T, calcNorm, spectrum, spatial_llh, log_interp_aeff, energy_llh, events)
    function signal_llh(params)
        Nex = params.Nex
        gamma = params.gamma
        E = params.E
        llh = zeros(Real, events.N)
        norm = calcNorm(Nex, gamma)
        for i = 1:events.N
            llh[i] += log(spectrum((norm = norm, E0=1e5, gamma=gamma, E=E[i], log10E=log10(E[i]), log10E0=log10(1e5))))
            llh[i] += spatial_llh[i]
            llh[i] += log_interp_aeff(log10(E[i]), sinDecPS)
            llh[i] += energy_llh(events.energy[i], log10(E[i]))
        end
        return llh
    end
    return signal_llh
end

function build_bg_llh(events, bg_norm, bg_llh)
    function return_func(params)
        Nex_bg = params.Nex_bg
        E = params.E
        log_Nex_bg = log(Nex_bg)
        llh = zeros(Real, events.N)
        for i = 1:events.N
            llh[i] = bg_llh[i] + bg_norm - log(E[i]) + log_Nex_bg
        end
        return llh
    end
    return return_func
end
```

```julia
energy_llh(3, ForwardDiff.Dual(3))
```

```julia
signal_llh = build_signal_llh(T, calcNorm, ps_coll.spectrum, spatial_llh, log_interp_aeff, energy_llh, events)
background_llh = build_bg_llh(events, bg_norm, bg_llh)
```

```julia
# signal_llh((Nex=ForwardDiff.Dual(1), E=fill(ForwardDiff.Dual(1e3), events.N), Nex_bg=ForwardDiff.Dual(10.), gamma=ForwardDiff.Dual(2.0)))
```

```julia
# background_llh((Nex=ForwardDiff.Dual(1), E=fill(ForwardDiff.Dual(1e3), events.N), Nex_bg=ForwardDiff.Dual(10.), gamma=ForwardDiff.Dual(2.0)))
```

```julia
function llh(params)
    event_llh = zeros(Real, events.N, 2)
    event_llh[:, 1] = signal_llh(params)
    event_llh[:, 2] = background_llh(params)
    llh = 0.
    llh = - params.Nex - params.Nex_bg
    for i=1:events.N
        llh += logsumexp(event_llh[i])
    end
    # return event_llh
    return llh
end
```

```julia
params = (Nex=10., Nex_bg=26., E=fill(1e3, events.N), gamma=2.0)
```

```julia
llh(params)
```

```julia
"""
function loglike(params)
    signal = signal_llh(params)
    # bg = background_llh(params)

    llh = -params.Nex # - params.Nex_bg
    # for i=1:length(signal)
    idxs = eachindex(signal)
    for i in idxs
        # llh += logsumexp(signal[i], bg[i])
        llh += signal[i]
    end
    llh
end
loglike(params)
"""
```

```julia
likelihood = let signal_llh = signal_llh, background_llh = background_llh
    logfuncdensity(function (params)
        signal = signal_llh(params)
        bg = background_llh(params)
        #llh = - params.Nex_bg + sum(bg)
        llh = -params.Nex - params.Nex_bg
        for i in eachindex(signal)
            llh += logsumexp(signal[i], bg[i])
        end
        return llh
    end)
end
```

```julia
prior = distprod(
    Nex = truncated(Normal(10, 5), 0., 100,),
    Nex_bg = truncated(Normal(25, 50), 0., 100.),
    gamma = truncated(Normal(2.3, 1.5), 1.2, 3.8),
    E = fill(Uniform(3e2, 3e8), events.N),
)
```

```julia
posterior = PosteriorMeasure(likelihood, prior)
```

```julia
context = set_batcontext(ad = ForwardDiff)
```

```julia
samples = bat_sample(
    posterior,
    TransformedMCMC(
        #proposal=HamiltonianMC(),
        proposal=RandomWalk(),
        nsteps = 10^4,
        nchains = 2,
        burnin = MCMCMultiCycleBurnin(max_ncycles=30, nsteps_per_cycle=1e5, nsteps_final=1e4),
    )
).result
```

```julia
using LazyReports
lazyreport(samples)
```

```julia
plot(
    samples, :(gamma), mean=true, std=true,
)
```

```julia
plot(
    samples, :(Nex), mean=true, std=true,
)
```

```julia
plot(
    samples, :(Nex_bg), mean=true, std=true,
)
```

```julia
copied
```

```julia
energies = flatview(samples.v.E[:, 1])
```

```julia
plot(log10.(energies[2, :]), mean=true, nbins=50)
```

```julia

```
