function R = GDFM_Extract(X, varargin)
%GDFM_EXTRACT Decompose multivariate data with four GDFM representations.
%
%   R = GDFM_Extract(X,'Method','TwoSided','q',2)
%   R = GDFM_Extract(X,'Method','OneSided','q',2,'r',2,'VARLags',2)
%   R = GDFM_Extract(X,'Method','Lee','q',2,'r',2,'VARLags',2)
%   R = GDFM_Extract(X,'Method','Frequency','q',2) % CPSD + SPCA
%   R = GDFM_Extract(X,'Method','Frequency','q',2, ...
%                    'StaticFactors',4) % optional reduced static rank
%
%   Input X is T-by-n (T observations, n observed series). Columns must be
%   finite, nonconstant, and approximately covariance stationary after any
%   user-specified trend/seasonality removal. The default standardizes each
%   column internally; disable using 'Standardize',false (centering remains
%   enabled). No forecasting or scenario generation occurs in this file.
%
%   Name/value options (case-insensitive):
%     Method        'TwoSided' (default), 'OneSided', 'Lee', 'Frequency'
%     q             dynamic eigenvector count; default min(2,n-1)
%     r             static factor count for OneSided/Lee; default q
%     m             covariance-lag half-width; default floor(sqrt(T))
%     f             frequency-grid half-width; default m (may differ)
%     CovMethod     'Bartlett_triang' (default), or 'Normal'
%     VARLags       finite VAR lag order p for OneSided/Lee (default 1)
%     Standardize   true (default) or false; always removes column means
%     StaticFactors  Number of static PCA directions (1..n, default n).
%                     Applies only to Method='Frequency'. A full-rank
%                     projection recovers the observed signal; a reduced
%                     rank projection retains only selected static modes.
%
%   Shared outputs: method, input, standardized, mean, scale, spectral,
%     factors, loadings, filters, identification, diagnostics. The
%     TwoSided/OneSided/Lee branches additionally populate common,
%     idiosyncratic, commonWithMean, and reconstruction (chi+xi).
%     Frequency does NOT populate those dynamic common-path fields.
%     Its separate frequencyTimeSeries structure contains the zero-lag
%     SPCA projection, scores, selected eigenvectors, and residual.
%
%   TWO-SIDED (Forni et al. 2000):
%     common = K(L) X; K(L) is a two-sided frequency projection filter.
%     factors is T-by-q obtained with the truncated extraction filter
%     B(L). loadings = A(w) is n-by-q-by-(2*f+1). The factor/loadings
%     decomposition is rotation/phase dependent, and the finite lag
%     factors*loadings filtered reconstruction does not in general equal
%     the direct K(L) estimate due to double filtering and truncation.
%
%   ONE-SIDED (Forni et al. 2005 - contemporaneous in-sample estimate):
%     factors=X*P; common=factors*loadings', where P uses generalized
%     eigenvectors of spectral common/idiosyncratic covariance matrices.
%     The generalized eigenproblem uses a regularized diagonal estimate
%     of the idiosyncratic covariance.
%     An optional modeling extension fits finite OLS VAR(p) factor dynamics.
%     The published forecasting step is NOT included here.
%
%   LEE (Lee & Baldick, IEEE TPWRS 2017, Eqs. 2-5, 13-16):
%     Estimates finite dynamic A(L), B(L) from spectral DPCA;
%     factors=B(L)X and common=A(L)factors. Estimates a finite VAR(p)
%     on the extracted dynamic factors using Yule-Walker covariance block
%     equations, and identifies epsilon=Z*delta using Cholesky.
%     Dynamic extraction includes past AND future observations (offline).
%     Innovations epsilon_t=C(L)f_t and whitened dynamic shocks delta_t
%     satisfy epsilon_t=Z*delta_t (Z is the Cholesky shock loading).
%     In the Lee branch, R.innovations stores epsilon_t, while R.shocks,
%     R.dynamicShocks, and R.whiteNoise all store delta_t. In OneSided,
%     only the unwhitened VAR innovations are calculated, so R.shocks
%     and R.dynamicShocks remain empty. The first p values are NaN.
%     Phase/rotation, truncation, factor scaling, VAR order, and Cholesky
%     ordering are finite-sample identification choices, not unique truths.
%     No scenario synthesis is performed.
%
%   FREQUENCY (Cho et al. 2026, Eqs. 8-15):
%     The core GDFM operation is S_X=S_chi+S_xi, with frequency-domain
%     A(w)=V_q(w), B(w)=A(w)^H, and f(w)=B(w)*S_X(w). f(w) is a
%     q-by-n-by-F spectral matrix, NOT a time-domain factor trajectory.
%
%     Default downstream static projection (separate from spectral GDFM):
%     IFFT of the unchanged S_X(w) -> estimated covariance Phi_X(k),
%     select zero lag Phi_X(0) -> static eigendecomposition P*D*P'.
%     The standardized time series is projected onto the retained r_s
%     static directions: X_SPCA = X_standardized*P_r*P_r'.
%     The static count r_s = StaticFactors is INDEPENDENT of q. By
%     default r_s=n, yielding X_SPCA=X to machine precision. For
%     r_s<n, X_SPCA is only a reduced-rank static projection.
%
%     Return this result as R.frequencyTimeSeries; do not confuse its
%     projected time series with a dynamically estimated common path.
%     The core frequency-domain CPSDs and filters never depend on r_s.
%     GDFM_FrequencySPCA(Snew,R.spectral.omega,R.standardized,r_s)
%     can also postprocess a modified CPSD; it does not synthesize
%     arbitrary new realizations matching that spectrum. At r_s=n,
%     the projection remains X even when the eigenvalues change.
%
%   For Lee, q == r is required for a q-dimensional dynamic factor path.
%   For methods with time-domain decomposition, common+idiosyncratic
%   equals standardized data BY CONSTRUCTION, not as an accuracy test. Means are restored
%   only in R.commonWithMean. See GDFM_RunAllExamples and GDFM_TestVariants.
%
%   References: Forni et al. (2000, 2005); Lee et al. (referenced as [18]
%   in Cho et al.); Cho, Zhu, Im, Lee, Baldick, IEEE TPWRS 2026, Table I,
%   Sections II.A-D. Dependencies: MATLAB base functionality only.

if ~isnumeric(X) || ~isreal(X) || ~ismatrix(X) || isempty(X) || ...
        any(~isfinite(X(:)))
    error('GDFM:Extract:Input','X must be a nonempty finite real matrix.');
end
[T,n] = size(X);
if n < 3 || T < 12
    error('GDFM:Extract:Size','Require at least 12 observations and 3 series.');
end
parser = inputParser;
parser.FunctionName = 'GDFM_Extract';
addParameter(parser,'Method','TwoSided',@(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(parser,'q',min(2,n-1),@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x==fix(x) && x>=1 && x<n);
addParameter(parser,'r',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x==fix(x) && x>=1 && x<n));
addParameter(parser,'m',floor(sqrt(T)),@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x==fix(x) && x>=1 && x<T);
addParameter(parser,'f',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x==fix(x) && x>=1));
addParameter(parser,'CovMethod','Bartlett_triang',@(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(parser,'VARLags',1,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x==fix(x) && x>=1);
addParameter(parser,'Standardize',true,@(x) islogical(x) && isscalar(x));
addParameter(parser,'StaticFactors',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x==fix(x) && x>=1 && x<=n));
parse(parser,varargin{:});
op = parser.Results;
if isempty(op.r), op.r = op.q; end
if isempty(op.f), op.f = op.m; end
if isempty(op.StaticFactors), op.StaticFactors = n; end
method = lower(char(op.Method));
switch method
    case {'twosided','two-sided'}, method = 'TwoSided';
    case {'onesided','one-sided'}, method = 'OneSided';
    case 'lee', method = 'Lee';
    case {'frequency','frequencydomain','frequency-domain'}, method = 'Frequency';
    otherwise
        error('GDFM:Extract:Method', ...
            'Method must be TwoSided, OneSided, Lee, or Frequency.');
end
if ~strcmp(method,'Frequency') && ...
        ~ismember('StaticFactors',parser.UsingDefaults)
    error('GDFM:Extract:FrequencyOnly', ...
        'StaticFactors applies only to Method=Frequency.');
end
% Avoid Fourier-grid aliasing at lag zero in the SPCA postprocessor.
% If the spectral grid is shorter than the covariance-lag support, the
% discrete inverse transform folds nonzero lag covariances onto lag zero.
if strcmp(method,'Frequency') && (2*op.f+1) <= op.m
    error('GDFM:Extract:FrequencyGrid', ...
        'Require 2*f+1 > m for unaliased zero-lag covariance in SPCA.');
end
if strcmp(method,'Lee') && op.r ~= op.q
    error('GDFM:Extract:LeeFactors', ...
        'Lee extraction requires r == q; dynamic factor count is q.');
end
if strcmp(method,'Lee') || strcmp(method,'OneSided')
    if T-op.VARLags <= op.r*op.VARLags
        error('GDFM:Extract:VARLength', ...
            'VARLags is too large relative to the sample size and r.');
    end
end

mu = mean(X,1);
scale = ones(1,n);
if any(std(X,0,1) <= 0)
    error('GDFM:Extract:Constant','X contains a constant series.');
end
if op.Standardize
    scale = std(X,0,1);
end
Z = bsxfun(@rdivide,bsxfun(@minus,X,mu),scale);
S = GDFM_SpectralCore(Z,op.q,op.m,op.f,char(op.CovMethod));

R = struct();
R.method = method;
R.input = X;
R.standardized = Z;
R.mean = mu;
R.scale = scale;
R.spectral = S;
R.factors = [];
R.loadings = [];
R.filters = struct();
R.factorSpectrum = [];
R.shocks = [];                % Identified dynamic shocks delta_t (Lee only)
R.dynamicShocks = [];         % Alias for shocks, delta_t (Lee only)
R.whiteNoise = [];            % Whitened VAR innovations delta_t (Lee only)
R.shockLoading = [];
R.noiseLoading = [];
R.innovations = [];
R.var = [];
R.identification = '';
R.common = []; R.idiosyncratic = []; R.commonWithMean = [];
R.idiosyncraticOriginal = []; R.reconstruction = [];
R.frequencyTimeSeries = [];  % Post-processing: separate from GDFM spectral factors

R.settings = struct('q',op.q,'r',op.r,'m',op.m,'f',op.f, ...
    'VARLags',op.VARLags,'CovMethod',char(op.CovMethod), ...
    'Standardize',op.Standardize, ...
    'StaticFactors',op.StaticFactors);

nf = 2*op.f+1;
A = S.eigvec;              % n-by-q-by-nf
B = complex(zeros(op.q,n,nf));
Fspec = complex(zeros(op.q,n,nf));
Cspec = complex(zeros(n,n,nf));
Kspec = complex(zeros(n,n,nf));
for j = 1:nf
    B(:,:,j) = A(:,:,j)';
    Fspec(:,:,j) = B(:,:,j)*S.full(:,:,j);
    Cspec(:,:,j) = A(:,:,j)*Fspec(:,:,j);
    Kspec(:,:,j) = A(:,:,j)*B(:,:,j);
end
R.factorSpectrum = Fspec;
R.filters.frequencyExtraction = B;
R.filters.frequencyLoading = A;
R.filters.frequencyProjector = Kspec;
R.filters.lags = -op.m:op.m;

switch method
    case 'TwoSided'
        % Finite two-sided time-domain estimator K(L)*Z.
        K = spectral_to_lag(Kspec,S.omega,R.filters.lags);
        R.filters.common = K;
        R.common = real(GDFM_ApplyFilter(Z,K,R.filters.lags));
        Bimp = spectral_to_lag(B,S.omega,R.filters.lags);
        Aimp = spectral_to_lag(A,S.omega,R.filters.lags);
        R.factors = real(GDFM_ApplyFilter(Z,Bimp,R.filters.lags));
        R.loadings = A;
        R.filters.factorExtraction = Bimp;
        R.filters.factorLoading = Aimp;
        R.identification = ['Spectral eigenvector phase/rotation is not ', ...
             'unique. Time-domain factors use a truncated two-sided ', ...
             'filter; composing truncated loading and extraction ', ...
             'filters need not reproduce K(L) exactly.'];
    case 'Frequency'
        % GDFM-specific operation: spectral eigenprojection only.
        R.factors = Fspec;     % q-by-n-by-F cross-spectral factor matrix
        R.loadings = A;       % n-by-q-by-F
        R.filters.dynamic = B; % q-by-n-by-F
        R.spectral.common = Cspec;
        R.spectral.idio = S.full-Cspec;
        R.identification = ['The Frequency branch provides S_common, ', ...
             'S_idio, A(w), B(w), and B(w)*S_X(w), not time-domain ', ...
             'factor trajectories. Optional time-series operations ', ...
             'are separate from spectral GDFM extraction.'];
        % Postprocessing is ALWAYS zero-lag SPCA (Cho et al., Eq. 14-15).
        % S.full has not been modified: it is the original estimated CPSD.
        % This is a static projection, not dynamic GDFM common extraction.
        H = GDFM_FrequencySPCA(S.full,S.omega,Z,op.StaticFactors);
        H.method = 'ZeroLagSPCA';
        H.projectedOriginal = bsxfun(@plus,bsxfun(@times, ...
            H.projectedStandardized,scale),mu);
        H.residualOriginal = bsxfun(@times,H.residualStandardized,scale);
        H.reconstructedOriginal = H.projectedOriginal+H.residualOriginal;
        R.frequencyTimeSeries = H;
    case 'OneSided'
        O = GDFM_OneSidedFactors(Z,S,op.r);
        R.common = O.common;
        R.factors = O.factors;         % T-by-r contemporaneous scores
        R.loadings = O.loadings;       % n-by-r static A0
        R.filters.staticExtraction = O.factorWeights;
        R.filters.staticLoading = O.loadings;
        R.static = O;
        V = GDFM_FactorVAR(O.factors,op.VARLags,false,'OLS');
        R.var = V;
        R.innovations = V.innovations;
        % No structural/whitened shocks are identified in this branch.
        R.noiseLoading = V.noiseLoading;
        R.identification = ['One-sided generalized static-factor estimator ', ...
              'with regularized diagonal idiosyncratic weights. ', ...
              'Finite OLS VAR dynamics are an extra assumption, ', ...
              'not the full Forni forecasting construction.'];
    case 'Lee'
        % Lee & Baldick (2017), Eqs. (2), (3), (13): A(L)B(L)X.
        L = GDFM_LeeDynamic(Z,S,op.q,R.filters.lags);
        R.common = L.common;
        R.factors = L.factors;          % T-by-q dynamic factor paths
        R.loadings = L.loadings;        % n-by-q-by-(2*m+1) dynamic A(L)
        R.filters.factorExtraction = L.extraction;
        R.filters.factorLoading = L.loadings;
        R.filters.frequencyExtraction = L.spectralFilter;
        R.filters.frequencyLoading = L.spectralLoading;
        R.filters.dynamic = L.spectralFilter;
        R.spectral.common = L.spectrumCommon;
        R.spectral.idio = S.full-L.spectrumCommon;
        R.lee = L;
        % Lee & Baldick (2017), Eqs. (4)-(5), (14)-(16): finite Yule-Walker
        % factor VAR followed by lower-triangular Cholesky identification.
        V = GDFM_FactorVAR(L.factors,op.VARLags,true,'YuleWalker');
        R.var = V;
        R.innovations = V.innovations;  % epsilon_t
        R.dynamicShocks = V.whitenedShocks;  % delta_t
        R.shocks = V.whitenedShocks;     % delta_t
        R.noiseLoading = V.noiseLoading;
        R.whiteNoise = V.whitenedShocks;
        R.shockLoading = V.shockLoading;
        R.identification = ['Lee-Baldick finite dynamic-filter model: ', ...
          'f=B(L)X, chi=A(L)f, Yule-Walker C(L), epsilon=Z delta. ', ...
          'Eigenvector phase, finite lags, VAR order and Cholesky ', ...
          'ordering are explicit identification/approximation choices.'];
    otherwise
        error('GDFM:Extract:Internal','Unexpected method.');
end

% Dynamic time-domain bookkeeping applies only when a common path was
% actually estimated. Static PCA is stored separately above.
R.diagnostics = struct();
R.diagnostics.commonFitScore = NaN;
R.diagnostics.maxReconstructionError = NaN;
if ~isempty(R.common)
    R.idiosyncratic = Z-R.common;
    R.commonWithMean = bsxfun(@plus,bsxfun(@times,R.common,scale),mu);
    R.idiosyncraticOriginal = bsxfun(@times,R.idiosyncratic,scale);
    R.reconstruction = R.commonWithMean+R.idiosyncraticOriginal;
    R.diagnostics.commonFitScore = 1-sum(R.idiosyncratic(:).^2)/sum(Z(:).^2);
    R.diagnostics.maxReconstructionError = max(abs(R.reconstruction(:)-X(:)));
end
R.diagnostics.maxSpectralIdentityError = max(abs(Cspec(:)-S.common(:)));
if strcmp(method,'Frequency')
    R.diagnostics.staticProjectionError = ...
        R.frequencyTimeSeries.relativeProjectionError;
    R.diagnostics.maxStaticDecompositionError = max(abs( ...
        R.frequencyTimeSeries.reconstructedOriginal(:)-X(:)));
end
if ~isempty(R.var)
    R.diagnostics.varSpectralRadius = R.var.spectralRadius;
    R.diagnostics.innovationStart = R.var.innovationStart;
    R.diagnostics.varEstimation = R.var.estimator;
    R.diagnostics.yuleWalkerMomentError = R.var.yuleWalkerMomentError;
    if strcmp(method,'Lee')
        R.diagnostics.leePhaseProjectionError = R.lee.phaseInvarianceError;
    end
    if strcmp(method,'Lee')
        valid = op.VARLags+1:T;
        delta = R.whiteNoise(valid,:);
        epsRec = (R.shockLoading*delta')';
        observedEps = R.innovations(valid,:);
        R.diagnostics.maxShockReconstructionError = ...
            max(abs(epsRec(:)-observedEps(:)));
    end
end
end

function H = spectral_to_lag(Hspec, omega, lags)
% Inverse discrete Fourier integral H_l = (1/nf) sum_w H(w) exp(i*l*w).
% No circular time convolution: GDFM_ApplyFilter uses zero padding.
[nOut,nIn,nf] = size(Hspec);
H = zeros(nOut,nIn,numel(lags));
for k = 1:numel(lags)
    temp = complex(zeros(nOut,nIn));
    for j = 1:nf
        temp = temp + Hspec(:,:,j)*exp(1i*lags(k)*omega(j));
    end
    H(:,:,k) = real(temp/nf);
end
end
