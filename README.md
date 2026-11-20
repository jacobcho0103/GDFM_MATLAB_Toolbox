# Generalized Dynamic Factor Model (GDFM) — MATLAB Toolbox

**Version 2.0.0**  |  **Author: Young-ho Cho**  |  **License: BSD 3-Clause**

A MATLAB implementation of generalized dynamic factor analysis for multivariate time series. The toolbox estimates common and idiosyncratic structure, spectral density matrices, factor coordinates, dynamic filters, and associated loading matrices. Four complementary formulations are provided, with reproducible examples of sensitivity analysis and illustrative scenario generation.

The primary MATLAB interface is **`GDFM_Extract`**. The tutorial and examples are documented in [examples/README.md](examples/README.md).

## Quick start

Requires MATLAB with base numerical and plotting functions. No external dataset or additional toolbox is needed for the included examples. Extract this folder, navigate to it in MATLAB, and run:

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

Input **`X`** is a real `T × n` matrix: rows are observations, columns are series. By default, columns are centered and standardized. For time-domain methods, the estimated decomposition is `Xstd = common + idiosyncratic`.

| Method | Principal calculation | Primary outputs |
|---|---|---|
| **TwoSided** — Forni et al. (2000) | `chi = K(L) Xstd` using a two-sided spectral projector | `common`, `idiosyncratic`, `filters.common`, finite factor coordinates |
| **OneSided** — Forni et al. (2005), contemporaneous estimate | `f = Xstd P`, `chi = f A0'` | Static `factors`, `loadings`, generalized eigenvalues; finite OLS-VAR option |
| **Lee** — Lee and Baldick (2017) | `f = B(L) Xstd`, `chi = A(L) f` | Dynamic factors and filters, Yule–Walker VAR, innovations and Cholesky-identified shocks |
| **Frequency** — Cho et al. (2026) | `SX = Scommon + Sidio` and `Scommon = A(w) B(w) SX` | Spectral matrices, `A(w)`, `B(w)`, matrix-valued spectral `f(w)` |

**Frequency-domain interpretation.** The `Frequency` method reports a `q × n × F` *spectral factor matrix*, rather than time-indexed factor trajectories. Its `R.common` field is empty by design. The separately calculated `R.frequencyTimeSeries` is a **zero-lag static PCA projection**, not a dynamically extracted GDFM common signal. That projection follows the covariance-based procedure in Cho et al. (2026), Eqs. (14)–(15).

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

The `Frequency` method requires `2*f+1 > m` to avoid zero-lag aliasing in its discrete inverse transform. For applications with known diurnal or seasonal patterns, remove deterministic trends and periodic components before fitting the stochastic GDFM.

## Illustrative scenario generation

Each method's example has a second part showing how an estimated component can be modified and how alternative time-series realizations can be constructed:

- **TwoSided:** vary the retained dynamic subspace and bootstrap residual blocks around the fixed common path.
- **OneSided:** intervene on a static factor, then simulate a fitted finite VAR and resample residual blocks.
- **Lee:** apply an identified shock impulse; generate new shock-driven factor trajectories through the VAR and lagged loadings.
- **Frequency:** modify a valid common CPSD, compare zero-lag SPCA projections, then separately synthesize stationary Gaussian time series from the modified spectrum.

These examples demonstrate stochastic constructions; their statistical properties depend on the fitted data and sampling assumptions. In particular, static PCA projection alone does not produce a new time-series realization matching a modified CPSD. The Frequency example therefore includes a distinct Fourier-domain Gaussian sampler. This package does not implement GAN training or a forecasting system.

## Files

```text
GDFM_MATLAB_Toolbox_v2.0.0/
├── GDFM_Extract.m                  Main interface
├── GDFM_SpectralCore.m             Spectral estimation and eigenvectors
├── GDFM_ApplyFilter.m              Finite matrix-valued filter
├── GDFM_OneSidedFactors.m          Contemporaneous static factors
├── GDFM_LeeDynamic.m               Dynamic factors and loadings
├── GDFM_FactorVAR.m                Factor VAR and shock identification
├── GDFM_FrequencySPCA.m            Separate zero-lag SPCA utility
├── GDFM_CheckInstallation.m        Installation checks
├── examples/                        Four methods and simulation helpers
├── tests/                           MATLAB numerical tests
├── LICENSE.txt
├── CITATION.cff
└── README.md
```

See [examples/README.md](examples/README.md) for the equations, figures, and interventions in each guided example.

## References and license

- Forni, Hallin, Lippi, and Reichlin, *The Generalized Dynamic-Factor Model: Identification and Estimation*, 2000.
- Forni, Hallin, Lippi, and Reichlin, *The Generalized Dynamic Factor Model: One-Sided Estimation and Forecasting*, 2005.
- Lee and Baldick, *Load and Wind Power Scenario Generation Through the Generalized Dynamic Factor Model*, IEEE Transactions on Power Systems, 2017. DOI: `10.1109/TPWRS.2016.2562718`.
- Cho, Zhu, Im, Lee, and Baldick, *Wind Power Scenario Generation Based on the Generalized Dynamic Factor Model and Generative Adversarial Network*, IEEE Transactions on Power Systems, 2026. DOI: `10.1109/TPWRS.2025.3615610`.

The implementation uses finite sample spectra and filters, with method-specific factor identification conventions. This distinction is documented directly in the function headers. The software is released under the [BSD 3-Clause License](LICENSE.txt). Citation metadata is provided in [CITATION.cff](CITATION.cff).
