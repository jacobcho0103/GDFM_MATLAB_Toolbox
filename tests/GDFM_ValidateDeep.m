function report = GDFM_ValidateDeep()
%GDFM_VALIDATEDEEP Independent numerical diagnostics for GDFM_Extract.
%   report = GDFM_ValidateDeep() uses simulated multivariate data with a
%   KNOWN common component and known Gaussian idiosyncratic noise.
%   Requires the GDFM_Extract.m dependencies in the same MATLAB path.
%
%   This complements GDFM_TestVariants. It checks whether estimated common
%   components are closer to ground truth than the noisy observations,
%   checks spectral/factor/filter algebra, and diagnoses recovered shocks.
%   It does NOT establish equivalence to the published Forni/Lee estimators;
%   that requires direct equation-by-equation and reference comparisons.
%   Only base MATLAB operations are used.
%
%   Important interpretation:
%     * Frequency returns a spectral GDFM decomposition and a separate
%       zero-lag SPCA projection, NOT a time-domain GDFM common path.
%     * OneSided uses generalized static factors, while Lee uses dynamic
%       two-sided B(L) factor extraction and A(L) common reconstruction.
%     * Lee shock covariance whitening is true algebraically because of
%       Cholesky; temporal whiteness is a separate empirical question.
%     * Factors and loadings are rotation/sign/phase non-unique, so their
%       raw values cannot be compared to simulated factors directly.
%
%   Reference signal: chi(t,:) = f(t,:)*A0' + f(t-1,:)*A1'.
%   Observation:      X(t,:)   = chi(t,:) + independent Gaussian noise.

methods = {'TwoSided','OneSided','Lee','Frequency'};
seeds = [37, 81, 124];
numSeeds = numel(seeds);
numMethods = numel(methods);
noiseBaseline = zeros(numSeeds,1);
relativeError = zeros(numSeeds,numMethods);
varianceFit = zeros(numSeeds,numMethods);
identityMax = zeros(numSeeds,numMethods);
frequencyIdentity = zeros(numSeeds,1);
shockLagCorrelation = zeros(numSeeds,1);
shockCovarianceError = zeros(numSeeds,1);
commonPairDifferences = zeros(numSeeds,2);

fprintf('\nGDFM independent recovery and component diagnostics\n');
fprintf('===================================================\n');
for s = 1:numSeeds
    rng(seeds(s));
    T = 1000; burnIn = 200; n = 12; q = 2; lagOrder = 2;
    a = linspace(0,2*pi,n)';
    A0 = [linspace(0.7,1.4,n)', 0.9*cos(a)];
    A1 = [0.22*sin(a), 0.18*cos(a+0.4)];
    u = randn(T+burnIn,q);
    trueF = zeros(T+burnIn,q);
    for t = 2:T+burnIn
        trueF(t,1) = 0.75*trueF(t-1,1)+u(t,1);
        trueF(t,2) = -0.45*trueF(t-1,2)+u(t,2);
    end
    trueF = trueF(burnIn+1:end,:);
    delayedF = [zeros(1,q); trueF(1:end-1,:)];
    trueChi = trueF*A0' + delayedF*A1';
    X = trueChi + 0.35*randn(T,n);

    signalScale = norm(bsxfun(@minus,trueChi,mean(trueChi,1)),'fro');
    noiseBaseline(s) = norm(X-trueChi,'fro')/signalScale;
    R = cell(1,numMethods);

    fprintf('\nSeed %d: noisy-observation error = %.4f\n', ...
        seeds(s),noiseBaseline(s));
    for k = 1:numMethods
        R{k} = GDFM_Extract(X, 'Method', methods{k}, ...
            'q',q,'r',q,'m',18,'f',25,'VARLags',lagOrder);
        if strcmp(methods{k},'Frequency')
            % SPCA full-rank projection is a downstream reconstruction.
            H=R{k}.frequencyTimeSeries;
            assert(isempty(R{k}.common));
            assert(norm(H.projectedOriginal-X,'fro')<1e-8);
            assert(norm(H.covariance0-R{k}.standardized'* ...
                R{k}.standardized/T,'fro')<1e-8);
            % Its full-rank error vs. true chi equals the noise baseline:
            % it is NOT a GDFM common-component recovery metric.
            relativeError(s,k)=NaN;
            varianceFit(s,k)=NaN;
            identityMax(s,k)=R{k}.diagnostics.maxStaticDecompositionError;
            fprintf('  %-10s spectral GDFM; full-rank SPCA input error %.2e\n', ...
                methods{k},norm(H.projectedOriginal-X,'fro'));
        else
            estChi = R{k}.commonWithMean;
            relativeError(s,k) = norm(estChi-trueChi,'fro')/signalScale;
            varianceFit(s,k) = R{k}.diagnostics.commonFitScore;
            identityMax(s,k) = R{k}.diagnostics.maxReconstructionError;
            assert(all(isfinite(estChi(:))), 'Nonfinite common component.');
            assert(identityMax(s,k) < 1e-10, 'Reconstruction identity fails.');
            fprintf('  %-10s relative true-common error = %.4f; fit = %6.2f%%\n', ...
                methods{k}, relativeError(s,k), 100*varianceFit(s,k));
            if relativeError(s,k) >= noiseBaseline(s)
                warning('GDFM:ValidateDeep:NoImprovement', ...
                    '%s did not improve over observed X for seed %d.', ...
                    methods{k},seeds(s));
            end
        end
    end

    % OneSided and Lee implement different time-domain common estimators.
    % The Frequency path is not equivalent to the TwoSided K(L) estimate.
    commonPairDifferences(s,1) = NaN;
    commonPairDifferences(s,2) = norm(R{2}.common-R{3}.common,'fro');

    % Frequency-domain identities from Eq. (8)-(13) in Cho et al. (2026).
    Freq = R{4};
    maxErr = 0;
    nf = size(Freq.spectral.full,3);
    for j = 1:nf
        A = Freq.loadings(:,:,j);
        B = Freq.filters.dynamic(:,:,j);
        Sj = Freq.spectral.full(:,:,j);
        Cj = Freq.spectral.common(:,:,j);
        fj = Freq.factors(:,:,j);
        maxErr = max(maxErr,norm(A'*A-eye(q),'fro'));
        maxErr = max(maxErr,norm(B-A','fro'));
        maxErr = max(maxErr,norm(fj-B*Sj,'fro'));
        maxErr = max(maxErr,norm(Cj-A*fj,'fro'));
        maxErr = max(maxErr,norm(Sj-Cj-Freq.spectral.idio(:,:,j),'fro'));
        maxErr = max(maxErr,norm(Sj-Sj','fro'));
    end
    frequencyIdentity(s) = maxErr;
    assert(maxErr < 1e-8, 'Frequency-domain identity fails.');

    % Check the finite two-sided lag-filter transpose symmetry K(-l)=K(l)'.
    two = R{1};
    K = two.filters.common;
    maxLagSymmetryError = 0;
    for j = 1:size(K,3)
        maxLagSymmetryError = max(maxLagSymmetryError, ...
            norm(K(:,:,j)-K(:,:,end+1-j)','fro'));
    end
    assert(maxLagSymmetryError < 1e-8, ...
        'Two-sided projection-lag symmetry fails.');

    % Lee branch uses a distinct dynamic factor filter A(L), B(L),
    % Yule-Walker VAR and Cholesky identification.
    Lee = R{3};
    valid = lagOrder+1:T;
    eps = Lee.innovations(valid,:);
    delta = Lee.whiteNoise(valid,:);
    Z = Lee.shockLoading;
    epsReconstructed = (Z*delta')';
    assert(norm(eps-epsReconstructed,'fro')/max(1,norm(eps,'fro')) < 1e-10);
    shockCovarianceError(s) = norm((delta'*delta)/numel(valid)-eye(q),'fro');
    assert(shockCovarianceError(s) < 1e-8, ...
        'Cholesky innovation covariance whitening fails.');

    % Sample serial correlations: small values support, but do not prove,
    % white-innovation behavior. No hard pass threshold for a stochastic test.
    maxCorr = 0;
    for d = 1:q
        v = delta(:,d)-mean(delta(:,d));
        for lag = 1:10
            rho = (v(1+lag:end)'*v(1:end-lag))/(v'*v);
            maxCorr = max(maxCorr,abs(rho));
        end
    end
    shockLagCorrelation(s) = maxCorr;
    fprintf('  Frequency identities max = %.2e; Lee shock covariance error = %.2e\n', ...
        maxErr,shockCovarianceError(s));
    fprintf('  Lee shock maximum |sample autocorrelation| at lags 1:10 = %.3f\n',maxCorr);
    fprintf('  Approximate single-lag 95%% noise band: +/- %.3f (not multiple-test corrected)\n', ...
        1.96/sqrt(numel(valid)));
end

fprintf('\nRecovery compared to the known simulated common component:\n');
for k = 1:numMethods
    if strcmp(methods{k},'Frequency')
        fprintf('  Frequency : spectral model; full-rank static SPCA reconstructs X, not chi.\n');
    else
        numImproved = sum(relativeError(:,k)<noiseBaseline);
        fprintf('  %-10s: improved on %d/%d seeds; mean normalized error %.4f\n', ...
            methods{k},numImproved,numSeeds,mean(relativeError(:,k)));
    end
end

report = struct('methods',{methods},'seeds',seeds, ...
    'observedError',noiseBaseline, 'estimatedCommonError',relativeError, ...
    'inSampleFit',varianceFit,'reconstructionIdentityError',identityMax, ...
    'spectralIdentityError',frequencyIdentity, ...
    'LeeShockCovarianceError',shockCovarianceError, ...
    'LeeShockMaxAutocorrelation',shockLagCorrelation, ...
    'estimatorPairDifferences',commonPairDifferences);

fprintf('\nAll deterministic numerical-identity checks passed.\n');
fprintf('Synthetic common-component recovery results are reported above.\n');
end
