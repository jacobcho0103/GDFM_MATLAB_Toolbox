# Generalized Dynamic Factor Model (GDFM) — MATLAB Toolbox

**Author:** Young-ho Cho, Duehee Lee, and Ross Baldick

A MATLAB implementation of generalized dynamic factor analysis for multivariate time series. The toolbox estimates common and idiosyncratic structure, spectral density matrices, factor coordinates, dynamic filters, and associated loading matrices. Four complementary formulations are provided, with reproducible examples of sensitivity analysis and illustrative scenario generation.

The primary MATLAB interface is **`GDFM_Extract`**. Optional seasonal adjustment uses MATLAB's built-in STL method through `GDFM_PreprocessSeasonality`, independently of GDFM extraction. The tutorial and examples are documented in [examples/README.md](examples/README.md).

## Quick start

The four GDFM methods and their examples use base MATLAB and require no external datasets or additional toolboxes. **Optional STL seasonal preprocessing requires MATLAB R2021b or newer** (`trenddecomp`, included in base MATLAB). Extract this folder, navigate to it in MATLAB, and run:

```matlab
addpath(genpath(pwd));
GDFM_CheckInstallation();     % checks extraction, SPCA, and examples
R = GDFM_Example_TwoSided();  % one fully worked method
```

To open all four demonstrations and their sensitivity/scenario figures:

```matlab
results = GDFM_RunAllExamples();
```

Each example returns a result structure containing estimated components and example-specific outputs. For figures-free runs, use `GDFM_Example_Lee(false,false)` (similarly for the other three examples).

## What the four methods calculate

Let `X` be a real `T × n` data matrix, with observations in rows and measured series in columns. `GDFM_Extract` centers every column and standardizes it by default. For the time-domain methods, the estimated decomposition satisfies `Xstd = common + idiosyncratic`. The mathematical expressions below use column-vector observations `x_t`; MATLAB stores the observations in rows.

| Method | Principal calculation | Main outputs |
|---|---|---|
| **TwoSided** [1] | $\widehat{\chi}_t=K(L)x_t$, with a two-sided spectral projection filter | Common and idiosyncratic time series, `filters.common` ($K(L)$), finite-filter factor coordinates |
| **OneSided** [2], in-sample static-factor estimate | $F=X_{\mathrm{std}}P_r$, $\widehat{X}^{\chi}=FA_0^{\mathsf T}$ | Static `factors`, `loadings`, generalized eigenvalues, optional finite OLS-VAR |
| **Lee** [3] | $\widehat{f}_t=B(L)x_t$, $\widehat{\chi}_t=A(L)\widehat{f}_t$ | Dynamic factors and filters, Yule–Walker VAR, innovations, Cholesky-normalized shocks |
| **Frequency** [4] | $S^X(\omega)=S^\chi(\omega)+S^\xi(\omega)$; $S^\chi(\omega)=A(\omega)B(\omega)S^X(\omega)$ | CPSD matrices, spectral loadings $A(\omega)$, filters $B(\omega)$, matrix-valued $f(\omega)$ |

The TwoSided method extracts the common signal through a finite two-sided approximation of the dynamic spectral projector. The OneSided implementation is a **contemporaneous in-sample static-factor estimator** associated with [2], not the complete one-sided forecasting algorithm. Lee's innovations are normalized by a Cholesky factor; covariance whitening does not establish uniquely identified structural shocks or temporal independence.

**Frequency-domain interpretation.** At each frequency, the leading $q$ spectral eigenvectors $V_q(\omega)$ define $A(\omega)=V_q(\omega)$, $B(\omega)=V_q(\omega)^H$, and $f(\omega)=B(\omega)S^X(\omega)$. Accordingly, $S^\chi(\omega)=A(\omega)f(\omega)$. The returned `R.factors` is a `q × n × F` **matrix-valued spectral quantity**, not a time-indexed factor trajectory. The `Frequency` method therefore does not populate `R.common`. Its `R.frequencyTimeSeries` output is a separate zero-lag covariance/static-PCA projection based on [4, eqs. (14)–(15)], interpreted in dimensionally consistent matrix form. This static projection is not a dynamic GDFM common-component estimator.

## Analyze your own data

```matlab
% X must be T-by-n, with finite, nonconstant columns.
[q, screen] = GDFM_SelectQ(X,'QMax',5,'m',16,'f',40);
R = GDFM_Extract(X,'Method','Lee','q',q, ...
    'r',q,'m',16,'f',40,'VARLags',2);
chi = R.commonWithMean;          % estimated common signal, original units
xi  = R.idiosyncraticOriginal;   % estimated residual, original units
```

For the spectral formulation and static time-series projection:

```matlab
R = GDFM_Extract(X,'Method','Frequency','q',q,'m',16,'f',40);
SX = R.spectral.full;                  % n-by-n-by-F CPSD
Scommon = R.spectral.common;           % common CPSD
Sidio = R.spectral.idio;               % idiosyncratic CPSD
A = R.loadings;                        % n-by-q-by-F spectral loading
B = R.filters.dynamic;                 % q-by-n-by-F dynamic filter
Fspectral = R.factors;                 % q-by-n-by-F spectral matrix
Xprojected = R.frequencyTimeSeries.projectedOriginal;
```

`StaticFactors` selects the number of static components for the Frequency projection **independently** of dynamic `q`. The default is all `n` static components, which gives the original input to numerical precision. With fewer components, the time series is a rank-reduced static projection.

```matlab
R = GDFM_Extract(X,'Method','Frequency','q',q, ...
    'm',16,'f',40,'StaticFactors',3);
```

### Optional seasonal preprocessing (MATLAB STL)

Seasonal adjustment is a separate preparation step, similar to factor-number screening; it does not change the four GDFM estimators. MATLAB's built-in `trenddecomp` provides **STL (seasonal-trend decomposition using Loess)** [5]. Apply it independently to each series before factor screening or GDFM estimation:

```matlab
% Example: hourly input with daily (24) and weekly (168) periodicities.
[Xadjusted, prep] = GDFM_PreprocessSeasonality(X,[24 168], ...
    'RemoveTrend',true);
[q, screen] = GDFM_SelectQ(Xadjusted,'QMax',5,'m',16,'f',40);
R = GDFM_Extract(Xadjusted,'Method','Lee','q',q,'r',q, ...
    'm',16,'f',40,'VARLags',2);

% Exact in-sample restoration of an adjusted series:
Xrestored = GDFM_RestoreSeasonality(Xadjusted,prep);
```

The default (`'RemoveTrend',false`) removes **seasonality only**. With `true`, STL's estimated long-term trend is also removed. `prep.seasonal`, `prep.trend`, and `prep.seasonalByPeriod` are accessible for inspection. Periods are specified in numbers of equally spaced observations and at least three cycles of the longest period are required. STL implements an *additive* decomposition; a transformed input may be preferable when seasonal amplitude scales with the signal level.

To restore same-horizon scenario paths `Yadjusted` of size `T x n x S`, use:

```matlab
YoriginalScale = GDFM_RestoreSeasonality(Yadjusted,prep);
```

The adjusted scenarios must first be in the **same units as `Xadjusted`**, not in GDFM-standardized units. The restoration utility adds the estimated seasonal (and optionally trend) components to each scenario; it is not a model for extrapolating seasonal components beyond the fitted sample. For a worked example with known synthetic seasonality and Lee GDFM, run:

```matlab
E = GDFM_Example_Seasonality();
% Or include it with the four method demonstrations:
results = GDFM_RunAllExamples('IncludeSeasonality',true);
```

The STL example is separate, so all four method examples run on MATLAB releases without `trenddecomp`. See the [official MATLAB STL documentation](https://www.mathworks.com/help/matlab/ref/double.trenddecomp.html).

### Factor-number selection

`GDFM_SelectQ.m` plots the averaged dynamic eigenspectrum and the adjacent eigenvalue ratios. It suggests a value of `q` at the largest ratio and returns `screen.eigenvalues`, `screen.ratios`, and `screen.recommendedQ`. Factor-number screening is a heuristic; inspect the eigenspectrum and consider the data's physical structure. The static rank for `OneSided` (`r`) and the Frequency postprocessing rank (`StaticFactors`) are distinct settings.

### Other parameters

| Parameter | Default | Meaning |
|---|---|---|
| `Method` | `'TwoSided'` | `'TwoSided'`, `'OneSided'`, `'Lee'`, or `'Frequency'` |
| `q` | `min(2,n-1)` | Number of leading dynamic spectral eigenvectors |
| `r` | `q` | Static factor count for OneSided; Lee uses `r=q` |
| `m` | `floor(sqrt(T))` | Covariance/filter lag half-width |
| `f` | `m` | Frequency-grid half-width; `F=2*f+1` |
| `CovMethod` | `'Bartlett_triang'` | Bartlett covariance weighting or `'Normal'` |
| `VARLags` | `1` | Finite factor VAR order, for OneSided and Lee |
| `Standardize` | `true` | Center and scale each column; `false` centers only |
| `StaticFactors` | `n` | Static PCA rank used only by Frequency |

The `Frequency` method requires `2*f+1 > m` to avoid zero-lag aliasing in its discrete inverse transform. For applications with known diurnal or seasonal patterns, use the optional STL preprocessing above or another appropriate seasonal-adjustment method before fitting the stochastic GDFM.

## Illustrative scenario generation

Each method's example has a second part showing how an estimated component can be modified and how alternative time-series realizations can be constructed:

- **TwoSided:** vary the retained dynamic subspace and bootstrap residual blocks around the fixed common path.
- **OneSided:** intervene on a static factor, then simulate a fitted finite VAR and resample residual blocks.
- **Lee:** apply a Cholesky-normalized shock impulse; generate new shock-driven factor trajectories through the VAR and lagged loadings.
- **Frequency:** modify a valid common CPSD, compare zero-lag SPCA projections, then separately synthesize stationary Gaussian time series from the modified spectrum.

These examples demonstrate stochastic constructions; their statistical properties depend on the fitted data and sampling assumptions. In particular, static PCA projection alone does not produce a new time-series realization matching a modified CPSD. The Frequency example therefore includes a distinct Fourier-domain Gaussian sampler. This package does not implement GAN training or a forecasting system.

## Files

```text
GDFM_MATLAB_Toolbox/
├── GDFM_Extract.m                  Main interface
├── GDFM_SpectralCore.m             Spectral estimation and eigenvectors
├── GDFM_ApplyFilter.m              Finite matrix-valued filter
├── GDFM_OneSidedFactors.m          Contemporaneous static factors
├── GDFM_LeeDynamic.m               Dynamic factors and loadings
├── GDFM_FactorVAR.m                Factor VAR and shock normalization
├── GDFM_FrequencySPCA.m            Separate zero-lag SPCA utility
├── GDFM_PreprocessSeasonality.m    Optional STL preprocessing
├── GDFM_RestoreSeasonality.m       Same-horizon seasonal restoration
├── GDFM_CheckInstallation.m        Installation checks
├── examples/                        Four methods, STL example, and simulation helpers
├── tests/                           MATLAB numerical tests including optional STL
├── CITATION.cff
└── README.md
```

See [examples/README.md](examples/README.md) for the equations, figures, and interventions in each guided example.

## References

[1] M. Forni, M. Hallin, M. Lippi, and L. Reichlin, “The generalized dynamic-factor model: Identification and estimation,” *Rev. Econ. Stat.*, vol. 82, no. 4, pp. 540–554, Nov. 2000, doi: [10.1162/003465300559037](https://doi.org/10.1162/003465300559037).

[2] M. Forni, M. Hallin, M. Lippi, and L. Reichlin, “The generalized dynamic factor model: One-sided estimation and forecasting,” *J. Amer. Stat. Assoc.*, vol. 100, no. 471, pp. 830–840, Sep. 2005, doi: [10.1198/016214504000002050](https://doi.org/10.1198/016214504000002050).

[3] D. Lee and R. Baldick, “Load and wind power scenario generation through the generalized dynamic factor model,” *IEEE Trans. Power Syst.*, vol. 32, no. 1, pp. 400–410, Jan. 2017, doi: [10.1109/TPWRS.2016.2562718](https://doi.org/10.1109/TPWRS.2016.2562718).

[4] Y.-H. Cho, H. Zhu, J. Im, D. Lee, and R. Baldick, “Wind power scenario generation based on the generalized dynamic factor model and generative adversarial network,” *IEEE Trans. Power Syst.*, vol. 41, no. 2, pp. 1135–1147, Mar. 2026, doi: [10.1109/TPWRS.2025.3615610](https://doi.org/10.1109/TPWRS.2025.3615610).

[5] R. B. Cleveland, W. S. Cleveland, J. E. McRae, and I. Terpenning, “STL: A seasonal-trend decomposition procedure based on Loess,” *J. Official Stat.*, vol. 6, no. 1, pp. 3–73, 1990.

Citation metadata for this software is provided in [CITATION.cff](CITATION.cff).
