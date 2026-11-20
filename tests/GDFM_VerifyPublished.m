function report = GDFM_VerifyPublished()
%GDFM_VERIFYPUBLISHED Equation-level tests for GDFM extraction methods.
%   report = GDFM_VerifyPublished()
%
%   Requires base MATLAB and GDFM_Extract plus its helper .m files.
%   Verifies deterministic definitions and synthetic recovery, NOT full
%   statistical consistency/identifiability for arbitrary input data.
%
%   Papers and equations:
%   - Forni et al. [1]: dynamic
%     principal-component projection of estimated cross spectrum.
%   - Forni et al. [2]: covariance-informed generalized static
%     factor estimator (the forecasting stage is out of scope).
%   - Lee and Baldick [3, eqs. (2)-(5), (13)-(16)]:
%       f=B(L)X; chi=A(L)f; eps=C(L)f; eps=Z*delta.
%   - Cho et al. [4, eqs. (6)-(13)]: frequency-domain
%     factorization of the observed cross power spectral density.
%
%   Different samples, different q/r and lag choices, and real-data checks
%   remain essential. The Lee A/B extraction is a finite two-sided filter,
%   not a causal streaming estimator. Shock normalization uses a Cholesky factor.

fprintf('\nGDFM published-equation verification\n');
fprintf('====================================\n');
methods = {'TwoSided','OneSided','Lee','Frequency'};
seeds = [37,81,124];
recovery = zeros(numel(seeds),4);
noiseError = zeros(numel(seeds),1);
maxDynamicIdentity = 0;
maxSpectralIdentity = 0;
maxYuleWalkerIdentity = 0;
maxCholeskyIdentity = 0;
maxOneSidedIdentity = 0;
for s = 1:numel(seeds)
    rng(seeds(s));
    T = 1000; burn = 200; n = 12; q = 2;
    t = linspace(0,2*pi,n)';
    A0 = [linspace(0.7,1.4,n)',0.9*cos(t)];
    A1 = [0.22*sin(t),0.18*cos(t+0.4)];
    noise = randn(T+burn,q);
    f = zeros(T+burn,q);
    for k = 2:T+burn
        f(k,1) = 0.75*f(k-1,1)+noise(k,1);
        f(k,2) = -0.45*f(k-1,2)+noise(k,2);
    end
    f = f(burn+1:end,:);
    chiTrue = f*A0' + [zeros(1,q); f(1:end-1,:)]*A1';
    X = chiTrue+0.35*randn(T,n);
    denom = norm(bsxfun(@minus,chiTrue,mean(chiTrue,1)),'fro');
    noiseError(s) = norm(X-chiTrue,'fro')/denom;
    R = cell(1,numel(methods));
    for k = 1:numel(methods)
        R{k}=GDFM_Extract(X,'Method',methods{k},'q',q,'r',q, ...
                          'm',18,'f',25,'VARLags',2);
        E = R{k};
        if strcmp(methods{k},'Frequency')
            recovery(s,k)=NaN; % CPSD representation is NOT a common path.
            H=E.frequencyTimeSeries;
            assert(isempty(E.common) && isempty(E.commonWithMean));
            assert(norm(H.projectedOriginal-X,'fro')<1e-8);
            assert(norm(H.covariance0-E.standardized'*E.standardized/T,'fro')<1e-8);
            assert(E.diagnostics.maxStaticDecompositionError<1e-10);
        else
            recovery(s,k) = norm(E.commonWithMean-chiTrue,'fro')/denom;
            assert(all(isfinite(E.common(:))));
            assert(E.diagnostics.maxReconstructionError < 1e-10);
            assert(recovery(s,k) < noiseError(s), ...
                'Common estimate did not improve over noisy input.');
        end
        assert(E.diagnostics.maxSpectralIdentityError < 1e-8);
    end
    two=R{1}; one=R{2}; lee=R{3}; freq=R{4};
    % Frequency is a CPSD factorization + zero-lag SPCA postprocessing.
    % It does NOT reuse the TwoSided K(L) common estimator.
    assert(isempty(freq.common));
    assert(norm(freq.frequencyTimeSeries.projectedOriginal-X,'fro') < 1e-8);
    % ONE-SIDED: verify generalized factor scores and static reconstruction.
    errOne=max([norm(one.factors-one.standardized*one.static.factorWeights,'fro'), ...
                norm(one.common-one.factors*one.loadings','fro')]);
    % Check the generalized eigenvalue equation Gamma_chi P = W P Lambda.
    Wdiag=max(one.static.weightsBeforeFloor,one.static.weightFloor);
    C0=one.static.commonCovariance0;
    P=one.static.factorWeights;
    eigvals=one.static.genEigenvalues;
    eigResidual=norm(C0*P-diag(Wdiag)*P*diag(eigvals),'fro') / ...
                max(1,norm(C0*P,'fro'));
    errOne=max(errOne,eigResidual);
    maxOneSidedIdentity=max(maxOneSidedIdentity,errOne);
    assert(errOne < 1e-8);
    % LEE: dynamic B(L) extraction and A(L) loading (not the static branch).
    lags=lee.filters.lags;
    fFromFilter=GDFM_ApplyFilter(lee.standardized,lee.filters.factorExtraction,lags);
    cFromFilter=GDFM_ApplyFilter(fFromFilter,lee.filters.factorLoading,lags);
    errLee=max([norm(fFromFilter-lee.factors,'fro'), ...
                norm(cFromFilter-lee.common,'fro')]);
    % Spectral adjoints imply time-lag duality B_l=A_{-l}^T.
    for il=1:numel(lags)
        errLee=max(errLee,norm(lee.filters.factorExtraction(:,:,il) - ...
            lee.filters.factorLoading(:,:,end+1-il)','fro'));
    end
    maxDynamicIdentity=max(maxDynamicIdentity,errLee);
    assert(errLee < 1e-8);
    assert(isequal(size(lee.loadings),[n,q,numel(lags)]));
    assert(strcmpi(lee.var.estimator,'YuleWalker'));
    assert(lee.var.yuleWalkerMomentError < 1e-9);
    maxYuleWalkerIdentity=max(maxYuleWalkerIdentity, ...
                               lee.var.yuleWalkerMomentError);
    assert(lee.diagnostics.leePhaseProjectionError < 1e-8);
    % Lee innovation and white-shock identities on their defined rows.
    valid=lee.var.innovationStart:T;
    eps=lee.innovations(valid,:);
    delta=lee.whiteNoise(valid,:);
    Z=lee.shockLoading;
    epsRec=delta*Z';
    errShock=max([norm(eps-epsRec,'fro')/max(1,norm(eps,'fro')), ...
      norm(delta'*delta/numel(valid)-eye(q),'fro')]);
    maxCholeskyIdentity=max(maxCholeskyIdentity,errShock);
    assert(errShock < 1e-8);
    assert(norm(eps-(lee.factors(valid,:)-lee.var.fitted(valid,:)), ...
                'fro') < 1e-9);
    % Independence is NOT tested by Cholesky; examine serial correlations.
    [maxAuto,maxCrossAuto] = local_lag_correlations(delta,10);
    % FREQUENCY: independent checks of the documented matrix identities.
    errSpec=0;
    for j=1:size(freq.spectral.full,3)
        A=freq.loadings(:,:,j);
        B=freq.filters.dynamic(:,:,j);
        Sj=freq.spectral.full(:,:,j);
        fj=freq.factors(:,:,j);
        Sc=freq.spectral.common(:,:,j);
        errSpec=max(errSpec,norm(A'*A-eye(q),'fro'));
        errSpec=max(errSpec,norm(B-A','fro'));
        errSpec=max(errSpec,norm(fj-B*Sj,'fro'));
        errSpec=max(errSpec,norm(Sc-A*fj,'fro'));
        errSpec=max(errSpec,norm(Sj-Sc-freq.spectral.idio(:,:,j),'fro'));
    end
    assert(errSpec<1e-8);
    maxSpectralIdentity=max(maxSpectralIdentity,errSpec);
    fprintf('Seed %d: noisy %.4f | Two %.4f | One %.4f | Lee %.4f | Frequency SPCA full-rank verified\n', ...
         seeds(s),noiseError(s),recovery(s,1:3));
    fprintf('  Lee shock max scalar autocorr=%.3f, max cross-lag corr=%.3f\n', ...
         maxAuto,maxCrossAuto);
end

% Separate oracle test: recover KNOWN VAR(2) coefficients from known factors.
% This validates the Yule-Walker implementation independently of GDFM factor
% identification. Cholesky-normalized shocks depend on factor ordering; structural identification is not established.
rng(901);
T=6000; burn=300;
M1=[0.55,0.06;-0.03,0.4];
M2=[0.12,0;0,0.1];
Ztrue=[1,0;0.25,0.75];
muTrue=[0.3;-0.2];
c=(eye(2)-M1-M2)*muTrue;
u=randn(T+burn,2)*Ztrue';
f=zeros(T+burn,2);
for t=3:T+burn
    f(t,:)=(c+M1*f(t-1,:)'+M2*f(t-2,:)'+u(t,:)')';
end
f=f(burn+1:end,:);
V=GDFM_FactorVAR(f,2,true,'YuleWalker');
coefError=max([norm(V.coefficient(:,:,1)-M1,'fro'), ...
               norm(V.coefficient(:,:,2)-M2,'fro')]);
assert(coefError < 0.12, ...
    'Yule-Walker failed known-VAR coefficient recovery.');
assert(V.yuleWalkerMomentError < 1e-9);
assert(V.spectralRadius < 1);
assert(norm(V.intercept-c','fro') < 0.12);
[shockAuto,shockCross] = local_lag_correlations(V.whitenedShocks(3:end,:),10);
fprintf('Known VAR(2) coefficient error: %.4f\n',coefError);
fprintf('Known VAR(2) whitened shocks: max auto %.3f, max cross %.3f\n', ...
    shockAuto,shockCross);

% Exact simulation identities are not equivalent to statistical validation;
% the uncertainty of finite samples is still present.
report = struct('methods',{methods},'seeds',seeds, ...
    'noiseBaseline',noiseError,'commonRecoveryError',recovery, ...
    'dynamicIdentityMax',maxDynamicIdentity, ...
    'oneSidedIdentityMax',maxOneSidedIdentity, ...
    'spectralIdentityMax',maxSpectralIdentity, ...
    'yuleWalkerIdentityMax',maxYuleWalkerIdentity, ...
    'choleskyIdentityMax',maxCholeskyIdentity, ...
    'knownVARCoefficientError',coefError, ...
    'knownVARShockMaxAuto',shockAuto, ...
    'knownVARShockMaxCrossLag',shockCross);
fprintf('PASS: All published-equation and known-VAR checks.\n');
fprintf('Equation-identity and synthetic VAR benchmarks completed.\n');
end

function [maxAuto,maxCrossLag] = local_lag_correlations(U,maxLag)
% Report maximum normalized lagged cross-correlations (diagnostic only).
[T,r]=size(U);
Uc=bsxfun(@minus,U,mean(U,1));
scale=sqrt(sum(Uc.^2,1));
maxAuto=0;
maxCrossLag=0;
for lag=1:maxLag
    for i=1:r
        for j=1:r
            rho=(Uc(lag+1:T,i)'*Uc(1:T-lag,j)) / ...
                 max(eps,scale(i)*scale(j));
            if i==j
                maxAuto=max(maxAuto,abs(rho));
            else
                maxCrossLag=max(maxCrossLag,abs(rho));
            end
        end
    end
end
end
