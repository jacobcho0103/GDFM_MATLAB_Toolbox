function GDFM_TestVariants()
%GDFM_TESTVARIANTS Smoke/regression tests for the extraction-only framework.
%   In MATLAB, add this folder to the path, then run GDFM_TestVariants.
%   Every test uses base MATLAB only. No test changes original GDFM.m.

rng(21);
T=340; n=7; r=2;
F=filter(1,[1,-0.5],randn(T,r));
A=[linspace(0.7,1.4,n)',cos((1:n)')];
X=F*A'+0.3*randn(T,n);
methods={'TwoSided','OneSided','Lee','Frequency'};
for i=1:numel(methods)
    R=GDFM_Extract(X,'Method',methods{i},'q',r,'r',r, ...
       'm',12,'f',17,'VARLags',2);
    if strcmp(R.method,'Frequency')
        assert(isempty(R.common) && isempty(R.idiosyncratic));
        assert(isempty(R.commonWithMean) && isempty(R.reconstruction));
        assert(isequal(size(R.frequencyTimeSeries.projectedOriginal),[T,n]));
        assert(R.diagnostics.maxStaticDecompositionError<1e-10);
    else
        assert(isequal(size(R.common),[T,n]));
        assert(isequal(size(R.idiosyncratic),[T,n]));
        assert(all(isfinite(R.common(:))));
        assert(R.diagnostics.maxReconstructionError<1e-10);
    end
    assert(R.diagnostics.maxSpectralIdentityError<1e-8);
    S=R.spectral;
    assert(isequal(size(S.common),[n,n,35]));
    assert(norm(S.common(:,:,18)-S.common(:,:,18)','fro')<1e-8);
    assert(norm(S.full(:)-S.common(:)-S.idio(:))<1e-8);
    switch R.method
        case 'TwoSided'
            assert(isequal(size(R.factors),[T,r]));
            assert(isequal(size(R.filters.factorExtraction),[r,n,25]));
            assert(isequal(size(R.loadings),[n,r,35]));
        case 'OneSided'
            assert(isequal(size(R.factors),[T,r]));
            assert(isequal(size(R.loadings),[n,r]));
            assert(isequal(size(R.noiseLoading),[r,r,3]));
            assert(isempty(R.dynamicShocks));
            assert(isempty(R.shocks));
            assert(isempty(R.whiteNoise));
            assert(all(all(isnan(R.innovations(1:2,:)))));
            assert(norm(R.common-R.factors*R.loadings','fro')<1e-8);
        case 'Lee'
            assert(isequal(size(R.factors),[T,r]));
            assert(isequal(size(R.loadings),[n,r,25]));
            assert(isequal(size(R.filters.factorExtraction),[r,n,25]));
            assert(strcmpi(R.var.estimator,'YuleWalker'));
            assert(R.var.yuleWalkerMomentError<1e-9);
            assert(norm(R.common-GDFM_ApplyFilter(R.factors,R.loadings,R.filters.lags),'fro')<1e-8);
            assert(isequal(size(R.shocks),[T,r]));
            assert(isequal(size(R.whiteNoise),[T,r]));
            assert(isequal(size(R.shockLoading),[r,r]));
            assert(all(all(isnan(R.whiteNoise(1:2,:)))));
            assert(all(all(isfinite(R.whiteNoise(3:end,:)))));
            assert(isequaln(R.shocks,R.whiteNoise));
            assert(isequaln(R.dynamicShocks,R.shocks));
            assert(norm(R.innovations(3:end,:)-R.shocks(3:end,:)*R.shockLoading','fro')<1e-8);
            assert(R.diagnostics.maxShockReconstructionError<1e-8);
        case 'Frequency'
            assert(isequal(size(R.factors),[r,n,35]));
            H=R.frequencyTimeSeries;
            assert(strcmp(H.method,'ZeroLagSPCA'));
            assert(H.staticFactors==n);
            assert(norm(H.covariance0-R.standardized'*R.standardized/T,'fro')<1e-8);
            assert(norm(H.projectedStandardized-R.standardized,'fro')<1e-8);
            assert(norm(H.projectedOriginal-X,'fro')<1e-8);
            assert(norm(H.projectedStandardized+H.residualStandardized-R.standardized,'fro')<1e-8);
            assert(norm(H.projection-eye(n),'fro')<1e-8);
            for k=[1,18,35]
                A0=R.loadings(:,:,k); f0=R.factors(:,:,k);
                assert(norm(A0*f0-R.spectral.common(:,:,k),'fro')<1e-8);
            end
            % Static rank is independent of dynamic q: reduced projection.
            Rsmall=GDFM_Extract(X,'Method','Frequency','q',r,'m',12, ...
                'f',17,'StaticFactors',r);
            Hr=Rsmall.frequencyTimeSeries;
            assert(Hr.staticFactors==r);
            assert(isequal(size(Hr.staticScores),[T,r]));
            assert(norm(Hr.projectedStandardized- ...
                Rsmall.standardized*Hr.retainedEigenvectors* ...
                Hr.retainedEigenvectors','fro')<1e-8);
            assert(norm(Hr.projectedStandardized-Rsmall.standardized,'fro')>1e-5);
            assert(norm(Rsmall.spectral.common(:)-R.spectral.common(:))<1e-8);
            assert(norm(Hr.reconstructedOriginal-X,'fro')<1e-8);
    end
    fprintf('PASS: %s\n',R.method);
end

% Chosen frequency grid f must not be confused with covariance lag m.
% The Frequency branch never substitutes TwoSided K(L) for its SPCA path.
R1=GDFM_Extract(X,'q',2,'m',8,'f',21);
assert(isequal(size(R1.spectral.eigvec),[n,2,43]));
assert(isequal(size(R1.filters.common),[n,n,17]));

% All invalid methods should fail clearly.
failed=false;
try
    GDFM_Extract(X,'Method','Unknown');
catch
    failed=true;
end
assert(failed);
fprintf('ALL GDFM extraction tests passed.\n');
end
