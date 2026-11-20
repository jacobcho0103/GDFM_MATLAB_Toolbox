function report = GDFM_StressTest()
%GDFM_STRESSTEST Independent oracle and robustness diagnostics for GDFM v3.
%   report = GDFM_StressTest();
%
%   Place alongside GDFM_Extract.m and its helper functions.
%   Base MATLAB only. No generation/scenario code is involved in the toolbox.
%
%   The TEST data are generated from a known two-factor dynamic model:
%       f_j(t) = a_j f_j(t-1) + u_j(t),  u_j ~ N(0,1)
%       chi(t) = A0*f(t) + A1*f(t-1)
%       X(t) = chi(t) + e(t),          e ~ N(0,sigma^2*I)
%
%   Verification dimensions:
%     1. Recovery of the known true common component, not chi+xi identity.
%     2. Frequency-domain CPSD compared with an INDEPENDENT analytic oracle.
%     3. Robustness over independent seeds, noise levels and sample lengths.
%     4. Yule-Walker shock whiteness as a diagnostic (not a strict pass test).
%     5. Finite-lag sensitivity of Lee and TwoSided common-component recovery.
%
%   IMPORTANT: Frequency-specific factors and dynamic eigenvectors are only
%   identifiable up to phase/rotation. Do not compare them entry-by-entry
%   with the generating factors/loadings. The oracle compares common spectra.
%   No finite test can prove estimator equivalence or asymptotic consistency.

methods = {'TwoSided','OneSided','Lee','Frequency'};
seeds = [11 23 47];
lengths = [800 1600];
noiseLevels = [0.25 0.50];
n = 10; q = 2; m = 18; fgrid = 25; varLags = 2;
a = [0.68; -0.38];
angles = linspace(0,2*pi,n)';
A0 = [linspace(0.8,1.25,n)', 0.75*cos(angles)];
A1 = [0.20*sin(angles), 0.14*cos(angles+0.3)];

nCases = numel(seeds)*numel(lengths)*numel(noiseLevels);
errors = NaN(nCases,numel(methods));
noisyErrors = NaN(nCases,1);
spectralErrors = NaN(nCases,1);
shockMaxAuto = NaN(nCases,1);
shockMaxCross = NaN(nCases,1);
shockAdjustedCutoff = NaN(nCases,1);
meta = NaN(nCases,3); % [seed, T, noise sd]
counter = 0;

fprintf('\nGDFM oracle and robustness verification\n');
fprintf('=======================================\n');
for iT = 1:numel(lengths)
    T = lengths(iT);
    for iN = 1:numel(noiseLevels)
        sigma = noiseLevels(iN);
        for iSeed = 1:numel(seeds)
            seed = seeds(iSeed);
            counter = counter+1;
            rng(seed+10000*iT+100*iN,'twister');
            [X,truth] = local_oracle_data(T,n,q,sigma,a,A0,A1);
            denominator = norm(truth-mean(truth,1),'fro');
            noisyErrors(counter) = norm(X-truth,'fro')/denominator;
            meta(counter,:) = [seed,T,sigma];
            results = cell(1,numel(methods));
            for k = 1:numel(methods)
                R = GDFM_Extract(X,'Method',methods{k},'q',q,'r',q, ...
                    'm',m,'f',fgrid,'VARLags',varLags);
                results{k} = R;
                if strcmp(methods{k},'Frequency')
                    % Frequency has only spectral GDFM outputs and a
                    % separate full-rank static PCA projection of X.
                    errors(counter,k)=NaN;
                    assert(isempty(R.common) && isempty(R.commonWithMean));
                    assert(norm(R.frequencyTimeSeries.projectedOriginal-X,'fro') < 1e-8);
                    assert(R.diagnostics.maxStaticDecompositionError < 1e-8);
                    assert(norm(R.frequencyTimeSeries.covariance0- ...
                        R.standardized'*R.standardized/T,'fro') < 1e-8);
                else
                    errors(counter,k) = norm(R.commonWithMean-truth,'fro')/denominator;
                    assert(isequal(size(R.common),size(X)),'Wrong common shape');
                    assert(all(isfinite(R.common(:))),'Nonfinite common estimate');
                    assert(R.diagnostics.maxReconstructionError < 1e-9, ...
                        'Reconstruction inconsistency');
                    assert(errors(counter,k) < noisyErrors(counter), ...
                        '%s did not improve on the noisy observation baseline',methods{k});
                end
            end
            % Oracle analytic common spectral density for each angular frequency.
            % Here spectra use sum of autocovariances, with no 1/(2*pi).
            S = results{4}.spectral;
            scale = results{4}.scale;
            oracleErrorSquared = 0;
            oracleNormSquared = 0;
            for j = 1:numel(S.omega)
                w = S.omega(j);
                Hw = A0+A1*exp(-1i*w);
                d = 1./abs(1-a*exp(-1i*w)).^2;
                Strue = Hw*diag(d)*Hw';
                Strue = bsxfun(@rdivide,bsxfun(@rdivide,Strue,scale'),scale);
                Se = S.common(:,:,j);
                oracleErrorSquared = oracleErrorSquared+norm(Se-Strue,'fro')^2;
                oracleNormSquared = oracleNormSquared+norm(Strue,'fro')^2;
                % The measured spectral matrices must obey these symmetries.
                assert(norm(S.full(:,:,j)-S.full(:,:,j)','fro')<1e-9, ...
                    'Non-Hermitian CPSD');
                jm = numel(S.omega)+1-j;
                assert(norm(S.full(:,:,jm)-conj(S.full(:,:,j)),'fro')<1e-9, ...
                    'CPSD conjugate-frequency symmetry broken');
            end
            spectralErrors(counter) = sqrt(oracleErrorSquared/oracleNormSquared);
            delta = results{3}.whiteNoise(varLags+1:end,:);
            [ma,mc] = local_shock_lag_correlation(delta,10);
            shockMaxAuto(counter)=ma;
            shockMaxCross(counter)=mc;
            % Approximate Bonferroni-adjusted two-sided normal cutoff for
            % q*q*10 cross/auto-lag checks. Diagnostic only, not a theorem.
            nComparisons = q*q*10;
            shockAdjustedCutoff(counter) = ...
                sqrt(2)*erfcinv(0.05/nComparisons)/sqrt(size(delta,1));
            fprintf('T=%4d sd=%.2f seed=%3d | noise %.3f | two %.3f | one %.3f | Lee %.3f | spectral oracle %.3f\n', ...
                T,sigma,seed,noisyErrors(counter),errors(counter,1:3),spectralErrors(counter));
        end
    end
end

% Sensitivity across admissible retained filter lag windows on a fixed data
% set. Lower error is not guaranteed to be monotonic as m varies.
rng(991,'twister');
[X,truth] = local_oracle_data(1200,n,q,0.35,a,A0,A1);
denominator = norm(truth-mean(truth,1),'fro');
lagOptions = [10 18 26];
lagSensitivity = zeros(numel(lagOptions),2);
for k=1:numel(lagOptions)
    rTwo = GDFM_Extract(X,'Method','TwoSided','q',q,'r',q, ...
        'm',lagOptions(k),'f',36);
    rLee = GDFM_Extract(X,'Method','Lee','q',q,'r',q, ...
        'm',lagOptions(k),'f',36,'VARLags',varLags);
    lagSensitivity(k,1) = norm(rTwo.commonWithMean-truth,'fro')/denominator;
    lagSensitivity(k,2) = norm(rLee.commonWithMean-truth,'fro')/denominator;
end
fprintf('\nLag-width sensitivity (same data, f=36):\n');
for k=1:numel(lagOptions)
    fprintf('m=%2d | two-sided %.4f | Lee %.4f\n', ...
        lagOptions(k),lagSensitivity(k,1),lagSensitivity(k,2));
end

fprintf('\nMean common recovery errors: two %.4f | one %.4f | Lee %.4f\n', ...
    mean(errors(:,1:3),1));
fprintf('Frequency is spectral; full-rank zero-lag SPCA is not a common-component estimator.\n');
fprintf('Mean analytic common-spectrum oracle error: %.4f\n',mean(spectralErrors));
fprintf('Lee shock max |auto| range: %.3f--%.3f\n', ...
    min(shockMaxAuto),max(shockMaxAuto));
fprintf('Approximate multiplicity-adjusted cutoff range: %.3f--%.3f\n', ...
    min(shockAdjustedCutoff),max(shockAdjustedCutoff));
fprintf('PASS: finite outputs, reconstruction identities, Hermitian/conjugate spectra, and recovery-versus-noise tests.\n');
fprintf('Spectral-oracle, filter sensitivity, and shock diagnostics are reported above.\n');

report = struct('methods',{methods},'caseMetadata',meta, ...
    'noisyError',noisyErrors,'commonRecoveryError',errors, ...
    'trueSpectralError',spectralErrors,'leeShockMaxAuto',shockMaxAuto, ...
    'leeShockMaxCrossLag',shockMaxCross, ...
    'leeShockAdjustedCutoff',shockAdjustedCutoff, ...
    'lagOptions',lagOptions,'lagSensitivityTwoAndLee',lagSensitivity);
end

function [X,chi] = local_oracle_data(T,n,q,sigma,a,A0,A1)
assert(q==2 && numel(a)==2 && size(A0,1)==n);
burn=300;
U=randn(T+burn,q);
F=zeros(T+burn,q);
for t=2:T+burn
    F(t,:) = F(t-1,:).*a'+U(t,:);
end
chiAll = F*A0'+[zeros(1,q);F(1:end-1,:)]*A1';
chi = chiAll(burn+1:end,:);
X = chi+sigma*randn(T,n);
end

function [maxAuto,maxCross] = local_shock_lag_correlation(U,maxLag)
[T,q]=size(U);
Uc=U-mean(U,1);
s=sqrt(sum(Uc.^2,1));
maxAuto=0; maxCross=0;
for lag=1:maxLag
    for i=1:q
        for j=1:q
            rho=(Uc(lag+1:T,i)'*Uc(1:T-lag,j))/max(eps,s(i)*s(j));
            if i==j
                maxAuto=max(maxAuto,abs(rho));
            else
                maxCross=max(maxCross,abs(rho));
            end
        end
    end
end
end
