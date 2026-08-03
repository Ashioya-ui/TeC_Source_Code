# Fork contributions

Two additions to `simonefatichi/TeC_Source_Code`, both benchmarked against
published data.

## 1. Making T&C FLEX-ready

FLEX/FLORIS launches September 2026 — 300 m SIF, 27-day repeat, in tandem with
Sentinel-3C. ESA's cal/val AO is open.

**T&C already models SIF.** `photosynthesis_biochemical.m` carries the Lee et al.
(2015, GCB) block: `Jfe` → `fiP` → `dls` → `kn` → `fiF` → `F755nm` in
W m⁻² sr⁻¹ µm⁻¹, which is the FLORIS retrieval band.
`Canopy_Resistence_An_Evolution.m` line 164 scales sunlit and shaded, and
`SIF_H`/`SIF_L` propagate to `MAIN_FRAME`.

**One step is missing.** Line 164 gives canopy *emission*. A satellite sees what
*escapes*.

`T&C_Code/SIF_Escape.m` adds it. Sunlit leaf density at depth `L` is
`exp(-Kopt·L)`, shaded the complement; escape toward zenith `θv` is `exp(-Kv·L)`:

```
fesc_sun = [(1-exp(-(Kopt+Kv)·LAI))/(Kopt+Kv)] / [(1-exp(-Kopt·LAI))/Kopt]
fesc_shd = [(1-exp(-Kv·LAI))/Kv - (1-exp(-(Kopt+Kv)·LAI))/(Kopt+Kv)]
           / [LAI - (1-exp(-Kopt·LAI))/Kopt]
```

These are the integrals already evaluated for the nitrogen profile in
`Canopy_Resistence_An_Evolution.m` lines 55–56, with viewing extinction `Kv` in
place of `Knit`. The denominators reduce to `LAI·Fsun` and `LAI·Fshd` **exactly**
as defined there — verified to 0.00e+00 over LAI 0.1–8.

`T&C_Code/SIF_Output.m` resamples onto FLORIS sampling using true local solar
time.

### Benchmark

Against Zhang, Joiner, Alemohammad, Zhou & Gentine (2018), *Biogeosciences* 15,
5779–5800 — CSIF evaluated against GPP at **40 FLUXNET tier-1 towers**:

| | modelled, Zurich, Vmax=55 | observed, 40 towers |
|---|---|---|
| slope (g C m⁻² d⁻¹ per mW m⁻² nm⁻¹ sr⁻¹) | **10.29** | 11.91 – 68.59 |
| r² | **0.896** | 0.01 – 0.93, median 0.64 |
| RMSE (g C m⁻² d⁻¹) | **1.86** | mean 1.67 |

The slope is 14% below the lower bound; r² sits at the top of the observed range
and RMSE is comparable. Nothing in the chain is fitted to this relationship, so
it is a prediction, not a calibration.

```
python Validation/benchmark_sif_gpp.py --lai 4 --vmax 55
```

### The result for cal/val

The slope is a strong function of `Vmax` and a weak function of LAI:

| Vmax | slope | r² | | LAI | slope | bulk fesc |
|---|---|---|---|---|---|---|
| 20 | 4.66 | 0.801 | | 1 | 7.37 | 0.915 |
| 40 | 8.25 | 0.866 | | 2 | 7.66 | 0.840 |
| 55 | 10.29 | 0.896 | | 4 | 8.08 | 0.712 |
| 80 | 13.19 | 0.930 | | 6 | 8.40 | 0.611 |
| 120 | 17.14 | 0.961 | | 8 | 8.63 | 0.529 |

LAI barely moves it because the escape fraction falls (0.915 → 0.529) roughly in
step with the GPP increase and the two largely cancel. `Vmax` moves it 3.7× over
a 6× range.

Two consequences. **SIF alone cannot constrain GPP without independent knowledge
of `Vmax`** — an apparent between-site slope difference may be a `Vmax`
difference. And **SIF and GPP jointly constrain `Vmax`**: at a tower measuring
both, the observed slope inverts to a `Vmax` estimate, with `r²` rising
monotonically alongside it. That is a usable cal/val target.

Note that `Vmax` alone cannot span the full observed 11.91–68.59: 3.7× over a
physiological `Vmax` range against 5.8× observed. Chlorophyll content, biome and
canopy structure carry the rest.

## 2. Forcing preparation from flux-tower data

T&C needs six radiation variables; a FLUXNET or ICOS tower gives total shortwave
only. `T&C_Code/Radiation_Partition.m` builds all six, and
`Forcing_Prep/prepare_forcing.py` produces a complete T&C `.mat` from tower CSV.

Two things found while calibrating against the shipped Zurich forcing:

**The timestamp convention is worth 19% of radiation RMSE and is undocumented.**
The forcing is stamped UTC while `DeltaGMT = 1`. Treating the stamp as local time
leaves 6,722 hours with Rsw > 20 W m⁻² and the sun below the horizon, and 3,771
hours with Rsw exceeding the extraterrestrial irradiance with the sun well up.
The correction *is* applied — by `t_bef = -0.67; t_aft = 1.67` in
`prova_Rural_Zurich.m`, whose window centres at +1.17 h, matching the physical
optimum of +1.10 h to 0.07 h. But those two constants carry no comment and are
site-specific, so copying the driver to a site with local-time stamps silently
imports Zurich's offset. `calibrate_hour_offset.m` determines it from the data.

**The visible fraction is a property of the stream, not the site.** In the shipped
bands, `SAB1/(SAB1+SAB2) = 0.389 ± 0.151` and `SAD1/(SAD1+SAD2) = 0.537 ± 0.101`
— 15 points apart, which is Rayleigh scattering. Any scheme applying one visible
fraction to both streams cannot reproduce the four bands: per-stream fractions
reconstruct them to RMSE 1.37 W m⁻², one shared fraction leaves 29.6 regardless.
Erbs + a fixed 0.45, and Weiss & Norman (1985), both make that assumption.

Out-of-sample (fit 1981–2004, tested 2005–2012), aggregate RMSE over the six
radiation variables falls **48.9%**.

## Limits

`Validation/chain.py` is a transcription of the MATLAB for benchmarking, not the
authoritative implementation. The escape formulation is single-scattering with a
two-stream albedo correction; a full treatment needs SCOPE-style radiative
transfer. The benchmark is one site, one PFT, constant LAI, against modelled
rather than measured SIF — a tower with a co-located spectrometer (DE-Hai) would
test it properly. The radiation coefficients are fitted at one mid-latitude
continental site and are an extrapolation elsewhere.
