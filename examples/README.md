# GDFM worked examples

Four method-specific MATLAB demonstrations illustrate the fitted quantities and downstream scenario constructions for the GDFM methods. Each of the four method examples uses the same reproducible stationary synthetic data so the decompositions can be compared consistently. A separate optional fifth example demonstrates STL seasonal adjustment before GDFM.

## Run

From the release root directory:

```matlab
addpath(genpath(pwd));
results = GDFM_RunAllExamples();
```

Or run one example:

```matlab
R = GDFM_Example_TwoSided();
R = GDFM_Example_OneSided();
R = GDFM_Example_Lee();
R = GDFM_Example_Frequency();

% Optional preprocessing demonstration, MATLAB R2021b or newer:
E = GDFM_Example_Seasonality();
```

To run the four method demonstrations followed by the optional seasonal example:

```matlab
results = GDFM_RunAllExamples('IncludeSeasonality',true);
```

A no-figure run, useful for automated tests, is `R = GDFM_Example_Frequency(false,false);`. The first logical argument controls the factor-screening figure; the second controls method and scenario figures. Each example returns the standard `GDFM_Extract` result structure, extended with `exampleSensitivity`, `exampleScenarios`, and `exampleScenarioStats`. The Frequency example also returns `exampleSpectralSynthesis`.

## Optional preparation: STL seasonality removal

`GDFM_PreprocessSeasonality.m` wraps MATLAB's established `trenddecomp(...,'stl',periods)` function based on STL [5]. Like `GDFM_SelectQ`, it is a stand-alone preparation utility, not part of the four GDFM estimators. Its inputs are an equally spaced `T x n` data matrix and one or more periods (measured in observations). The example uses a known 24-observation seasonal pattern and a known trend, then fits Lee's GDFM representation to the adjusted signal.

```matlab
% Adjust a multivariate signal before GDFM analysis:
[Xadj, prep] = GDFM_PreprocessSeasonality(X,24,'RemoveTrend',true);
R = GDFM_Extract(Xadj,'Method','TwoSided','q',2,'m',16,'f',40);

% Restore the seasonal and trend components after any method's
% generated scenarios have been converted to adjusted-input units:
Yraw = GDFM_RestoreSeasonality(Yadj,prep);  % T-by-n-by-S
```

`RemoveTrend=false` by default retains the estimated long-term trend. The preprocessing result stores the estimated trend, individual seasonal components, their sum, and the removed component. The restoration utility supports only the same time horizon as the fitted STL decomposition. It adds deterministic components to all scenarios; it does not forecast them. For an end-to-end worked example, use `GDFM_Example_Seasonality()` (requires MATLAB R2021b or newer). The four GDFM method examples remain unchanged and have no STL dependency.

## Shared data and factor screening

`GDFM_DemoData.m` creates 960 observations of 12 correlated series from **two known dynamic factors**, lagged loadings, and idiosyncratic noise. It returns the data and the generating common and residual components, enabling independent estimation-error checks. The examples use a fixed random seed.

`GDFM_SelectQ.m` estimates the dynamic eigenspectrum of the CPSD and plots both average eigenvalues and adjacent eigenvalue ratios. For the included synthetic data, the largest ratio suggests `q=2`, matching the data generator. The selected `q` determines the retained dynamic eigenspace. For `OneSided`, the static count `r` can differ from `q`; for Frequency, `StaticFactors` is a separate static projection rank.

### 1. TwoSided — dynamic projection and residual resampling

**Decomposition [1].** The common signal is estimated by a two-sided matrix-valued filter, `chi = K(L) Xstd`; the idiosyncratic signal is the difference `Xstd - chi`. The filter uses lagged and leading observations.

**First figure.** Compare the observed, true common, and estimated common signals; compare true and estimated residuals; examine extracted finite-filter factor coordinates and the lag-dependent Frobenius norm of `K(L)`.

**Sensitivity.** Re-estimate the spectral projector using `q+1` rather than `q`, and compare the common signal and `K(L)` coefficients. This changes the retained subspace while preserving the projector construction.

**Scenario example.** Retain the estimated common path and generate new residuals by jointly moving-block resampling the estimated idiosyncratic series. Add each resampled residual path to the common path, then restore the original mean and scale. This is a conditional-residual scenario construction, rather than a newly sampled common-factor path.

```matlab
R = GDFM_Example_TwoSided();
K = R.filters.common;
Kmodified = R.exampleSensitivity.filterAlternative;
Y = R.exampleScenarios;       % T-by-n-by-5
```

### 2. OneSided — contemporaneous factors and factor dynamics

**Decomposition [2].** The in-sample common signal is `chi = f A0'`, with contemporaneous factors `f = Xstd P` estimated using a generalized covariance eigenproblem. This is the static-factor reconstruction associated with the one-sided framework; it is distinct from a complete causal forecasting procedure.

**First figure.** Compare common and residual paths with synthetic truth; inspect static factors, loading matrix `A0`, selected generalized eigenvalues, and residuals of an additional finite VAR fit.

**Sensitivity.** Scale the first estimated static factor by 1.30 and hold the loading matrix fixed to show its effect across observed series.

**Scenario example.** Fit an OLS VAR to the factor paths, resample factor-innovation rows, propagate the VAR, reconstruct through `A0`, and add jointly resampled residual blocks. The finite VAR is an additional simulation assumption.

```matlab
R = GDFM_Example_OneSided();
f = R.factors;
A0 = R.loadings;
Y = R.exampleScenarios;
```

### 3. Lee — dynamic filters, innovations, and shock propagation

**Extraction [3].** Estimate finite two-sided filters `A(L)` and `B(L)` from the leading spectral eigenvectors. Extract dynamic factors through `f = B(L) Xstd` and common signals through `chi = A(L) f`. Estimate finite VAR parameters by Yule–Walker moment equations, then use a Cholesky factor to normalize the innovations as `epsilon_t = Z delta_t`.

**First figure.** Compare the common and residual paths with synthetic truth; inspect dynamic factor trajectories, the lagged loading/filter matrices, VAR innovations, and covariance-whitened shocks.

**Sensitivity.** Add a three-unit impulse to the first Cholesky-normalized shock and propagate its effect through the fitted VAR and dynamic loadings.

**Scenario example.** Generate independent Gaussian shocks, apply the fitted shock-loading matrix, simulate factors using the VAR, reconstruct with `A(L)`, and add block-bootstrap residuals. The shock representation is conditional on the fitted factor coordinates and Cholesky ordering.

```matlab
R = GDFM_Example_Lee();
f = R.factors;
AL = R.filters.factorLoading;
BL = R.filters.factorExtraction;
delta = R.shocks;
Y = R.exampleScenarios;
```

### 4. Frequency — spectral factors and separate time-series synthesis

**GDFM spectral decomposition [4].** Estimate the multivariate CPSD and separate it into the leading `q` spectral directions and the remaining directions:

```text
Sx(w) = Scommon(w) + Sidio(w)
A(w)  = Vq(w)
B(w)  = A(w)^H
f(w)  = B(w) Sx(w)
Scommon(w) = A(w) f(w)
```

The `f(w)` output is a `q-by-n-by-F` *spectral factor matrix*, rather than a time-domain factor path.

**Separate zero-lag SPCA.** Inverse Fourier transform the total CPSD, select lag-zero covariance, and use its leading static eigenvectors to project observed data: `Xstatic = Xstd P_r P_r'`. The default `StaticFactors=n` reproduces the input; a reduced static rank gives a different projection. This operation is downstream of GDFM spectral extraction and is not a dynamic GDFM common-component estimate.

**First figure.** Compare full and reduced SPCA projections and residuals; inspect zero-lag covariance eigenvalues, the total/common CPSD power, magnitudes of `A(w)`, `B(w)`, `f(w)`, and the covariance heatmap.

**Sensitivity.** Apply a smooth positive even frequency-dependent weighting to the common CPSD while keeping the idiosyncratic CPSD fixed. Compare original and modified spectra, covariance eigenvalues, and reduced-rank SPCA projections. Full-rank SPCA still projects to the original data because the complete static basis spans the observation space.

**Scenario example.** Use a separate stationary Gaussian Fourier sampler on the modified CPSD. It interpolates spectral matrices to FFT frequencies, samples conjugate-symmetric Gaussian Fourier coefficients, and transforms them to real time-series paths. Expected second-order statistics approximate the specified spectrum; finite paths exhibit sampling variation. This sampler is separate from SPCA and does not use a GAN.

```matlab
R = GDFM_Example_Frequency();
SX = R.spectral.full;
Scommon = R.spectral.common;
A = R.loadings;
B = R.filters.dynamic;
fSpectral = R.factors;
Xstatic = R.frequencyTimeSeries.projectedOriginal;
Smodified = R.exampleSensitivity.modifiedSpectrum;
Y = R.exampleScenarios;       % T-by-n-by-30
```

## Scenario outputs and diagnostic interpretation

All four examples provide `R.exampleScenarios` (`T × n × S`) in the input's original units. `R.exampleScenarioStats` includes relative cross-series covariance discrepancy and mean absolute lag-one correlation discrepancy for each scenario. These diagnostics illustrate specific statistical properties, not a single overall scenario-quality score. Domain applications may additionally require checks for periodicity, marginal and ramp-rate distributions, operating limits, extreme events, and dependence across longer time horizons.

The shared demonstration helpers are `GDFM_DemoResidualBootstrap.m`, `GDFM_DemoSimulateVAR.m`, `GDFM_DemoScenarioDiagnostics.m`, and `GDFM_DemoSpectralSynthesis.m`. Their inputs and outputs are described in MATLAB help comments.

## Automated checks

```matlab
addpath(genpath(pwd));
GDFM_CheckInstallation();
```

This checks all four extraction methods, frequency/SPCA reconstruction, and each sensitivity/scenario example against deterministic numerical properties. On MATLAB R2021b or newer, it also runs `GDFM_TestSeasonalPreprocessing` for the optional STL workflow. The `tests/` folder also contains deeper synthetic recovery, factor VAR, and spectral-oracle diagnostics.

Method references [1]–[4] and the STL reference [5] are listed in the [toolbox README](../README.md).
