function R = GDFM_Example_Frequency(showScreening,showPlots)
%GDFM_EXAMPLE_FREQUENCY Spectral GDFM plus separate zero-lag SPCA projection.
%
%   R = GDFM_Example_Frequency();
%   R = GDFM_Example_Frequency(false); % omit repeated q-screening figure
%   R = GDFM_Example_Frequency(false,false); % compute without figures
%
%   Frequency-domain GDFM [4, eqs. (8)-(13)]:
%       SX(w) = Scommon(w) + Sidio(w)
%       A(w) = Vq(w), B(w) = A(w)', f(w) = B(w)*SX(w)
%       Scommon(w) = A(w)*f(w)
%   Here f(w) is a q-by-n cross-spectral matrix, NOT a factor time series.
%
%   Downstream reconstruction (Eqs. 14-15; NOT the GDFM decomposition):
%       SX(w) -> IFFT -> PhiX(k) -> PhiX(0) -> eig(PhiX(0))
%       Xstatic = Xstd * P_r * P_r'
%   r = n by default, which recovers the input exactly; r < n gives a
%   reduced-rank static projection. Static r is independent of dynamic q.
%
%   This SPCA operation is not the TwoSided K(L)*X estimate and cannot
%   generate arbitrary new time series from a modified CPSD by itself.
%
%   Requires core GDFM_Extract.m, GDFM_FrequencySPCA.m, GDFM_DemoData.m,
%   and GDFM_SelectQ.m. No external data or special MATLAB toolbox.

%   Sensitivity example: apply a PSD-preserving frequency weighting.
%   Compare reduced-rank SPCA with separate stochastic Fourier synthesis.

if nargin < 1, showScreening = true; end
if nargin < 2, showPlots = true; end
if ~islogical(showPlots) || ~isscalar(showPlots)
    error('GDFM:Example:PlotFlag','showPlots must be a logical scalar.');
end
if ~islogical(showScreening) || ~isscalar(showScreening)
    error('GDFM:Example:PlotFlag','showScreening must be a logical scalar.');
end
exampleDir = fileparts(mfilename('fullpath'));
addpath(exampleDir);
addpath(fileparts(exampleDir));
[X,truth] = GDFM_DemoData('Seed',2026);
m=16; f=40;
[q,~] = GDFM_SelectQ(X,'QMax',5,'m',m,'f',f, ...
    'Plot',showScreening);

% Dynamic q controls the spectral decomposition. StaticFactors controls
% only the downstream SPCA projection; full rank is the default.
R = GDFM_Extract(X,'Method','Frequency','q',q,'m',m,'f',f);
H = R.frequencyTimeSeries;
[nSeries,~] = size(R.spectral.full(:,:,1));
staticReduced = min(q,nSeries-1);
Hreduced = GDFM_FrequencySPCA(R.spectral.full,R.spectral.omega, ...
    R.standardized,staticReduced);
Xreduced = bsxfun(@plus,bsxfun(@times, ...
    Hreduced.projectedStandardized,R.scale),R.mean);
R.exampleReducedSPCA = Hreduced;

omega=R.spectral.omega;
nf=numel(omega);
k0=f+1+round(f/2); % nonzero frequency: show complex factor magnitudes
A=R.loadings(:,:,k0);
B=R.filters.dynamic(:,:,k0);
Fspec=R.factors(:,:,k0);
maxFilterErr=0; maxFactorErr=0; maxCommonErr=0;
traceFull=zeros(1,nf); traceCommon=zeros(1,nf);
for k=1:nf
    Ak=R.loadings(:,:,k); Bk=R.filters.dynamic(:,:,k);
    Fk=R.factors(:,:,k);
    maxFilterErr=max(maxFilterErr,norm(Bk-Ak','fro'));
    maxFactorErr=max(maxFactorErr,norm(Fk-Bk*R.spectral.full(:,:,k),'fro'));
    maxCommonErr=max(maxCommonErr, ...
        norm(R.spectral.common(:,:,k)-Ak*Fk,'fro'));
    traceFull(k)=real(trace(R.spectral.full(:,:,k)));
    traceCommon(k)=real(trace(R.spectral.common(:,:,k)));
end

% The zero-lag covariance from the IFFT agrees with the unweighted
% sample covariance when the spectral grid avoids lag aliasing (F>m).
Csample = R.standardized'*R.standardized/size(X,1);
covarianceErr = norm(H.covariance0-Csample,'fro');
fullError = norm(H.projectedOriginal-X,'fro')/max(1,norm(X,'fro'));
reducedError = norm(Xreduced-X,'fro')/max(1,norm(X,'fro'));

fprintf('\n[Frequency] dynamic q=%d (synthetic generating q=%d)\n',q,truth.q);
fprintf('  Spectral outputs: A %s, B %s, f %s\n', ...
    mat2str(size(R.loadings)),mat2str(size(R.filters.dynamic)), ...
    mat2str(size(R.factors)));
fprintf('  Spectral checks (all %d frequencies): |B-A^H| %.2e, |f-B*S| %.2e, |Scommon-A*f| %.2e\n', ...
    nf,maxFilterErr,maxFactorErr,maxCommonErr);
fprintf('  Zero-lag IFFT covariance vs sample covariance: %.3e\n',covarianceErr);
fprintf('  SPCA projection: full static rank=%d, relative input error %.3e\n', ...
    H.staticFactors,fullError);
fprintf('  SPCA projection: reduced static rank=%d, relative input error %.3e\n', ...
    staticReduced,reducedError);
fprintf('  The SPCA projection is separate from spectral GDFM extraction.\n');

ix=1:min(250,size(X,1));
if showPlots
figure('Name','Example 4: Frequency-domain GDFM and zero-lag SPCA','Color','w');
subplot(2,4,1);
plot(ix,X(ix,1),'k','LineWidth',1.1); hold on;
plot(ix,H.projectedOriginal(ix,1),'b--','LineWidth',1.0);
plot(ix,Xreduced(ix,1),'r','LineWidth',1.0);
legend('Observed','Full-rank SPCA','Reduced SPCA','Location','best');
xlabel('Time'); title('Static time-series projection'); grid on;
subplot(2,4,2);
plot(ix,X(ix,1)-H.projectedOriginal(ix,1),'b','LineWidth',1); hold on;
plot(ix,X(ix,1)-Xreduced(ix,1),'r','LineWidth',1);
legend('Full-rank residual','Reduced-rank residual','Location','best');
xlabel('Time'); title('SPCA residual (not GDFM idio)'); grid on;
subplot(2,4,3);
plot(1:numel(H.eigenvalues),H.eigenvalues,'ko-'); hold on;
yl = ylim; plot([staticReduced staticReduced],yl,'r--');
xlabel('Static PC index'); title('Zero-lag covariance eigenvalues'); grid on;
subplot(2,4,4);
plot(omega,traceFull,'k-'); hold on;
plot(omega,traceCommon,'b-','LineWidth',1.1);
legend('Trace S_X','Trace S_{common}','Location','best');
xlabel('Angular frequency'); title('Spectral decomposition'); grid on;
subplot(2,4,5);
imagesc(abs(A)); colorbar;
xlabel('Dynamic factor'); ylabel('Series'); title('|A(\omega)| spectral loading');
subplot(2,4,6);
imagesc(abs(B)); colorbar;
xlabel('Series'); ylabel('Dynamic factor'); title('|B(\omega)| spectral filter');
subplot(2,4,7);
imagesc(abs(Fspec)); colorbar;
xlabel('Series'); ylabel('Dynamic factor'); title('|f(\omega)| spectral factor matrix');
subplot(2,4,8);
imagesc(H.covariance0); colorbar;
xlabel('Series'); ylabel('Series'); title('IFFT zero-lag covariance');
end % decomposition and projection figures

%% Spectral sensitivity, SPCA projection, and stochastic scenario synthesis
% Frequency-domain GDFM extraction is unchanged. Here, illustrate a valid
% alternative (PSD) spectrum by multiplying the common CPSD by a smooth,
% positive, EVEN frequency weight. Retain the estimated idiosyncratic CPSD.
% This is spectral sensitivity, NOT the GAN-synthesized dynamic filter from
% [4]; it does not change the fitted GDFM spectral factors.
frequencyGain=1+1.5*exp(-0.5*((abs(omega)-0.8)/0.35).^2);
Smodified=complex(zeros(size(R.spectral.full)));
traceModified=zeros(1,nf);
maxPSDAdjustment=0;
% Bartlett/windowed empirical spectra can have small negative eigenvalues
% because finite-sample covariance estimates are not necessarily PSD.
% Enforce PSD BEFORE sampling, then mirror the conjugate frequency partner.
% This projection is downstream scenario preparation, not GDFM extraction.
for j=1:f+1
    Sj=R.spectral.idio(:,:,j)+ ...
        frequencyGain(j)*R.spectral.common(:,:,j);
    Sj=(Sj+Sj')/2;
    if j==f+1, Sj=real(Sj); end % zero-frequency CPSD must be real
    [U,D]=eig(Sj);
    SjPSD=U*diag(max(real(diag(D)),0))*U';
    SjPSD=(SjPSD+SjPSD')/2;
    if j==f+1, SjPSD=real(SjPSD); end
    maxPSDAdjustment=max(maxPSDAdjustment, ...
        norm(SjPSD-Sj,'fro')/max(1,norm(Sj,'fro')));
    Smodified(:,:,j)=SjPSD;
    if j<=f
        Smodified(:,:,nf+1-j)=conj(SjPSD);
    end
end
for j=1:nf
    traceModified(j)=real(trace(Smodified(:,:,j)));
end
% The zero-lag SPCA is a separate projection. Its eigenvectors can change
% under the modified CPSD. At full rank it still returns the original X.
HmodifiedReduced=GDFM_FrequencySPCA(Smodified,omega,R.standardized,staticReduced);
HmodifiedFull=GDFM_FrequencySPCA(Smodified,omega,R.standardized,nSeries);
XmodifiedSPCA=bsxfun(@plus,bsxfun(@times, ...
    HmodifiedReduced.projectedStandardized,R.scale),R.mean);
fullRankChange=norm(HmodifiedFull.projectedStandardized-R.standardized,'fro') ...
    /max(1,norm(R.standardized,'fro'));
reducedChange=norm(HmodifiedReduced.projectedStandardized- ...
    Hreduced.projectedStandardized,'fro')/max(1,norm(R.standardized,'fro'));

% To obtain NEW random waveforms rather than an SPCA projection, separately
% synthesize Gaussian Fourier coefficients from the modified positive-
% semidefinite CPSD. Expected second-order statistics approximate the
% target, but a finite realization has no exact CPSD guarantee.
oldState=rng;
cleanup=onCleanup(@() rng(oldState)); %#ok<NASGU>
rng(4204,'twister');
nScenarios=30;
[Ystandardized,synthInfo]=GDFM_DemoSpectralSynthesis( ...
    Smodified,omega,size(X,1),nScenarios);
scenarios=zeros([size(X),nScenarios]);
for s=1:nScenarios
    scenarios(:,:,s)=bsxfun(@plus,bsxfun(@times, ...
        Ystandardized(:,:,s),R.scale),R.mean);
end
stats=GDFM_DemoScenarioDiagnostics(X,scenarios);
R.exampleSensitivity=struct('frequencyGain',frequencyGain, ...
    'modifiedSpectrum',Smodified,'modifiedSPCA',HmodifiedReduced, ...
    'maxPSDAdjustment',maxPSDAdjustment, ...
    'relativeReducedProjectionChange',reducedChange, ...
    'relativeFullRankProjectionError',fullRankChange);
R.exampleScenarios=scenarios;
R.exampleSpectralSynthesis=synthInfo;
R.exampleScenarioStats=stats;
fprintf('  [Frequency extension] weighted common CPSD; zero-lag SPCA projection change %.3g (reduced) / %.3g (full).\n', ...
    reducedChange,fullRankChange);
fprintf('  Spectral preparation: largest relative PSD-cone correction %.3g.\n',maxPSDAdjustment);
fprintf('  Spectral-Gaussian sampling: %d scenarios; ensemble covariance-to-target error %.3g.\n', ...
    nScenarios,synthInfo.relativeEnsembleCovarianceError);
fprintf('  New paths are NOT created by zero-lag SPCA. Frequency interpolation and finite-sample differences remain.\n');
if showPlots
    ix=1:min(250,size(X,1));
    figure('Name','Frequency extension: modified CPSD, SPCA, spectral sampling','Color','w');
    subplot(2,3,1);
    plot(omega,traceFull,'k','LineWidth',1.1); hold on;
    plot(omega,traceModified,'r--','LineWidth',1.1);
    legend('Original total CPSD','Modified total CPSD','Location','best');
    xlabel('Angular frequency'); title('Smooth frequency-selective amplification'); grid on;
    subplot(2,3,2);
    plot(Hreduced.eigenvalues,'k-o'); hold on;
    plot(HmodifiedReduced.eigenvalues,'r--x');
    legend('Original','Modified','Location','best');
    xlabel('Static eigenvalue index'); title('IFFT zero-lag covariance eigenvalues'); grid on;
    subplot(2,3,3);
    plot(ix,Xreduced(ix,1),'b','LineWidth',1.1); hold on;
    plot(ix,XmodifiedSPCA(ix,1),'r--','LineWidth',1.1);
    legend('Original reduced SPCA','Modified reduced SPCA','Location','best');
    xlabel('Time'); title('SPCA is a deterministic projection'); grid on;
    subplot(2,3,4);
    plot(ix,X(ix,1),'k','LineWidth',1.2); hold on;
    plot(ix,squeeze(scenarios(ix,1,1:5)),'LineWidth',0.8);
    xlabel('Time'); title('Independent Gaussian spectral realizations'); grid on;
    subplot(2,3,5);
    imagesc(synthInfo.targetCovariance); colorbar;
    xlabel('Series'); ylabel('Series'); title('Interpolated target covariance');
    subplot(2,3,6);
    imagesc(synthInfo.meanEmpiricalCovariance); colorbar;
    xlabel('Series'); ylabel('Series'); title('Ensemble empirical covariance');
end

end
