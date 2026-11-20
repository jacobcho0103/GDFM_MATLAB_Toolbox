function report = GDFM_TestFrequencySPCA()
%GDFM_TESTFREQUENCYSPCA Verify spectral GDFM vs zero-lag SPCA postprocessing.
%
%   report = GDFM_TestFrequencySPCA();
%
%   Verifies the explicit separation:
%   (1) Frequency GDFM returns S_X=S_common+S_idio, A(w), B(w), f(w).
%   (2) Its ONLY time-series postprocessing path is IFFT->Phi(0)->SPCA.
%   (3) By default StaticFactors=n reconstructs X (not latent chi).
%   (4) Reduced rank projects X onto the top static eigenvectors.
%   (5) Rotating an otherwise valid CPSD changes the reduced projection,
%       but not the full-rank projection. Eigenvalue-only shifts have no
%       effect on the static projection.
%
%   These are deterministic matrix tests, NOT evidence that a modified
%   spectrum by itself generates a new time series consistent with that CPSD.

rng(103,'twister');
T=480; n=6; q=2; m=10; fgrid=18;
factor=filter(1,[1 -0.6],randn(T,q));
loadings=[linspace(1,2,n)' cos((1:n)')];
X=factor*loadings' + 0.25*randn(T,n) + ones(T,1)*(1:n);
R=GDFM_Extract(X,'Method','Frequency','q',q,'m',m,'f',fgrid);
H=R.frequencyTimeSeries;
assert(isempty(R.common) && isempty(R.idiosyncratic));
assert(isempty(R.reconstruction));
assert(strcmp(H.method,'ZeroLagSPCA'));
assert(H.staticFactors==n);
assert(norm(H.covariance0-R.standardized'*R.standardized/T,'fro')<1e-8);
assert(norm(H.projectedOriginal-X,'fro')<1e-8);
assert(norm(H.projectedStandardized+H.residualStandardized- ...
    R.standardized,'fro')<1e-8);

% Frequency decomposition identities on ALL frequencies.
maxSpectral=0;
for k=1:size(R.spectral.full,3)
    A=R.loadings(:,:,k); B=R.filters.dynamic(:,:,k);
    SX=R.spectral.full(:,:,k); FS=R.factors(:,:,k);
    maxSpectral=max(maxSpectral,norm(B-A','fro'));
    maxSpectral=max(maxSpectral,norm(FS-B*SX,'fro'));
    maxSpectral=max(maxSpectral,norm(R.spectral.common(:,:,k)-A*FS,'fro'));
    maxSpectral=max(maxSpectral,norm(SX-R.spectral.common(:,:,k)- ...
        R.spectral.idio(:,:,k),'fro'));
end
assert(maxSpectral<1e-9);

% StaticFactors is independently selectable and does not change spectra.
Rr=GDFM_Extract(X,'Method','Frequency','q',q,'m',m, ...
    'f',fgrid,'StaticFactors',q);
Hr=Rr.frequencyTimeSeries;
assert(norm(Rr.spectral.full(:)-R.spectral.full(:))<1e-10);
assert(norm(Hr.projectedStandardized- ...
    R.standardized*Hr.retainedEigenvectors*Hr.retainedEigenvectors','fro')<1e-8);
assert(norm(Hr.projectedStandardized-R.standardized,'fro')>1e-4);
assert(norm(Hr.reconstructedOriginal-X,'fro')<1e-8);

% Postprocess modified CPSD without changing spectral estimator.
Q=eye(n); theta=0.6;
Q([1 3],[1 3])=[cos(theta),-sin(theta); sin(theta),cos(theta)];
Sm=R.spectral.full;
for k=1:size(Sm,3)
    Sm(:,:,k)=Q*Sm(:,:,k)*Q';
end
Hmod=GDFM_FrequencySPCA(Sm,R.spectral.omega,R.standardized,q);
HmodFull=GDFM_FrequencySPCA(Sm,R.spectral.omega,R.standardized,n);
rotatedDifference=norm(Hmod.projectedStandardized-Hr.projectedStandardized,'fro');
assert(rotatedDifference>1e-4);
assert(norm(HmodFull.projectedStandardized-R.standardized,'fro')<1e-8);

% Uniform spectral eigenvalue shift does not alter eigenvectors / projection.
Sshift=R.spectral.full;
for k=1:size(Sshift,3)
    Sshift(:,:,k)=Sshift(:,:,k)+0.4*eye(n);
end
Hshift=GDFM_FrequencySPCA(Sshift,R.spectral.omega,R.standardized,q);
assert(norm(Hshift.projectedStandardized-Hr.projectedStandardized,'fro')<1e-7);

report=struct('spectralIdentityError',maxSpectral, ...
    'fullRankError',norm(H.projectedOriginal-X,'fro'), ...
    'reducedRankError',norm(Hr.projectedOriginal-X,'fro'), ...
    'rotatedSpectrumReducedProjectionChange',rotatedDifference, ...
    'modifiedSpectrumFullRankError', ...
       norm(HmodFull.projectedStandardized-R.standardized,'fro'));
fprintf('PASS: Frequency spectral identities (%.2e)\n',maxSpectral);
fprintf('PASS: IFFT zero-lag covariance / default full-rank SPCA (%.2e)\n', ...
    report.fullRankError);
fprintf('PASS: Independent reduced static rank (input difference %.3e)\n', ...
    report.reducedRankError);
fprintf('PASS: Modified CPSD rotates reduced projection (change %.3e)\n', ...
    rotatedDifference);
fprintf('PASS: Full-rank projection unchanged by modified CPSD\n');
fprintf('PASS: Uniform eigenvalue shift preserves static projection\n');
fprintf('ALL Frequency/SPCA interface tests passed.\n');
end
