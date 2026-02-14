---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.1
  kernelspec:
    display_name: hnu
    language: python
    name: python3
---

```python
from astropy.coordinates import SkyCoord
import astropy.units as u
from hierarchical_nu.source.parameter import Parameter
from hierarchical_nu.simulation import Simulation
from hierarchical_nu.fit import StanFit
from hierarchical_nu.priors import Priors
from hierarchical_nu.source.source import Sources, PointSource, DetectorFrame
from hierarchical_nu.utils.lifetime import LifeTime
from hierarchical_nu.events import Events
from hierarchical_nu.fit import StanFit
from hierarchical_nu.priors import Priors, LogNormalPrior, NormalPrior, LuminosityPrior, IndexPrior, FluxPrior, DifferentialFluxPrior
from hierarchical_nu.utils.roi import CircularROI, FullSkyROI
from hierarchical_nu.detector.icecube import IC86_II, IC86_I
from hierarchical_nu.detector.input import mceq
from icecube_tools.utils.data import Uptime
from hierarchical_nu.source.cosmology import *
import numpy as np
import ligo.skymap.plot
```

```python
# define high-level parameters
Parameter.clear_registry()
src_index = Parameter(2.1, "src_index", fixed=False, par_range=(1, 4))
diff_index = Parameter(2.52, "diff_index", fixed=False, par_range=(1, 4))
L = Parameter(1e47 * (u.erg / u.s), "luminosity", fixed=True, 
              par_range=(0, 1E60) * (u.erg/u.s))
#diffuse_norm = Parameter(2.26e-13 /u.GeV/u.m**2/u.s, "diffuse_norm", fixed=True, 
                         #par_range=(1e-14, 1e-11)*(1/u.GeV/u.s/u.m**2))
z = 0.3365
Enorm = Parameter(1E5 * u.GeV, "Enorm", fixed=True)
Emin = Parameter(1E2 * u.GeV, "Emin", fixed=True)
Emax = Parameter(1E9 * u.GeV, "Emax", fixed=True)
Emin_src = Parameter(Emin.value, "Emin_src", fixed=True)
Emax_src = Parameter(Emax.value, "Emax_src", fixed=True)
#Emin_diff = Parameter(Emin.value, "Emin_diff", fixed=True)
#Emax_diff = Parameter(Emax.value, "Emax_diff", fixed=True)
```

```python
Emin_det = Parameter(1e1 * u.GeV, "Emin_det", fixed=True)
```

```python
# Single PS for testing and usual components
ra = np.deg2rad(77.35) * u.rad
dec = np.deg2rad(5.7) * u.rad
width = np.deg2rad(6) * u.rad
txs = SkyCoord(ra=ra, dec=dec, frame="icrs")
point_source = PointSource.make_powerlaw_source(
    "test", dec, ra, L, src_index, z, Emin_src, Emax_src, DetectorFrame,
)

my_sources = Sources()
my_sources.add(point_source)
my_sources.add_background(IC86_II)
```

```python
lt = LifeTime()
obs_time = lt.lifetime_from_dm(IC86_II)
```

```python
obs_time
```

```python
MJD_min=56917
MJD_max=57113

# roi = CircularROI(txs, 5 * u.deg, MJD_min=MJD_min, MJD_max=MJD_max, apply_roi=True)
roi = FullSkyROI()
```

```python
event_types = [IC86_II]
events = Events.from_ev_file(
    *event_types)
print(events.N)
```

```python
fit = StanFit(my_sources, IC86_II, events, obs_time)
```

```python
fit.precomputation()
```

```python
inputs = fit._get_fit_inputs()
```

```python
len(inputs["bg_llh"])
bg_llh = (np.array(inputs["bg_llh"]) + np.log(obs_time[IC86_II].to_value(u.s)))
np.savetxt("ic86_ii_bg_llh.csv", bg_llh)
```

```python
sim = Simulation(my_sources, event_types, lifetime)
```

```python
sim.precomputation()
print(sim._get_expected_Nnu(sim._get_sim_inputs()))
print(sim._expected_Nnu_per_comp)
```

```python
point_source.parameters
```

```python
E(1)
```

```python
hubble_factor(1.)
```

```python
comoving_distance(1.)
```

```python
luminosity_distance(1.)
```

```python
redshift(1 * u.Mpc)
```

```python

```
