function R = GDFM_Example_TwoSided(showScreening,showPlots)
%GDFM_EXAMPLE_TWOSIDED Example 1: Forni (2000) two-sided decomposition.
%
%   R = GDFM_Example_TwoSided();
%   R = GDFM_Example_TwoSided(false); % omit repeated q-screening figure
%   R = GDFM_Example_TwoSided(false,false); % compute without figures
%
%   Steps: simulate stationary multivariate data with known common signal;
%   screen the number q of dynamic eigenvectors; estimate chi and xi;
%   inspect the two-sided common filter K(L) and illustrative factor paths.
%   K(L) uses future AND past data, so it is not an online forecasting filter.
%   Needs GDFM_Extract.m and its five core helpers on the MATLAB path,
%   plus GDFM_DemoData.m and GDFM_SelectQ.m in this folder.

%   Scenario example: re-estimate K(L) with different q, then generate conditional
%   scenarios by jointly block-resampling the estimated residuals.

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
addpath(fileparts(exampleDir));  % Parent folder holds the core MATLAB files.
if exist('GDFM_Extract','file')~=2
    error('GDFM:Example:Path','Add the GDFM core folder to the MATLAB path.');
end
[X,truth] = GDFM_DemoData('Seed',2026);
m = 16; f = 40;
[q,screen] = GDFM_SelectQ(X,'QMax',5,'m',m,'f',f,'Plot',showScreening); %#ok<NASGU>
R = GDFM_Extract(X,'Method','TwoSided','q',q,'m',m,'f',f);

trueCommon = bsxfun(@plus,truth.common,truth.mean);
relativeError = norm(R.commonWithMean-trueCommon,'fro')/norm(truth.common,'fro');
fprintf('\n[TwoSided] q=%d (synthetic truth q=%d)\n',q,truth.q);
fprintf('  Estimated common shape: %d x %d\n',size(R.common));
fprintf('  Relative error vs known common: %.4f\n',relativeError);
fprintf('  Two-sided projector K(L): %d x %d x %d lags\n',size(R.filters.common));
fprintf('  Common + idiosyncratic identity: %.2e\n',R.diagnostics.maxReconstructionError);
fprintf('  Factor paths are rotation/phase-dependent finite-filter coordinates.\n');

ix = 1:min(250,size(X,1));
K = R.filters.common;
klag = zeros(1,size(K,3));
for j=1:size(K,3)
    klag(j) = norm(K(:,:,j),'fro');
end
if showPlots
figure('Name','Example 1: TwoSided - decomposition and filter','Color','w');
subplot(2,2,1);
plot(ix,X(ix,1),'Color',[0.65,0.65,0.65]); hold on;
plot(ix,trueCommon(ix,1),'k--','LineWidth',1.2);
plot(ix,R.commonWithMean(ix,1),'b','LineWidth',1.1);
legend('Observed','True common','Estimated common','Location','best');
xlabel('Time index'); ylabel('Series 1'); title('Two-sided common component'); grid on;
subplot(2,2,2);
plot(ix,truth.idiosyncratic(ix,1),'k--'); hold on;
plot(ix,R.idiosyncraticOriginal(ix,1),'r');
legend('True idiosyncratic','Estimated idiosyncratic','Location','best');
xlabel('Time index'); title('Idiosyncratic component'); grid on;
subplot(2,2,3);
plot(ix,R.factors(ix,1:min(2,q)),'LineWidth',1);
xlabel('Time index'); ylabel('Factor score');
title(sprintf('Extracted two-sided factor paths (q=%d)',q)); grid on;
subplot(2,2,4);
plot(R.filters.lags,klag,'o-','LineWidth',1.2);
xlabel('Filter lag (negative = future)'); ylabel('Frobenius norm');
title('Two-sided common filter K(L)'); grid on;
end % decomposition figures

%% Scenario and sensitivity demonstration: sensitivity of K(L), then conditional scenario generation
% Vary the dynamic subspace dimension while rebuilding K(L) through the
% estimator; arbitrary coefficient edits could destroy its projector form.
qAlternative = min(size(X,2)-1,q+1);
Ralt = GDFM_Extract(X,'Method','TwoSided','q',qAlternative,'m',m,'f',f);
commonDifference = norm(Ralt.common-R.common,'fro')/max(1,norm(R.common,'fro'));
Kalternative = Ralt.filters.common;
altLagNorm = zeros(1,numel(Ralt.filters.lags));
for j=1:numel(altLagNorm)
    altLagNorm(j)=norm(Kalternative(:,:,j),'fro');
end

% A fixed fitted common path plus JOINT moving-block resampled idiosyncratic
% deviations: conditional/residual-bootstrap scenarios, not new GDFM factors.
oldState=rng;
cleanup=onCleanup(@() rng(oldState)); %#ok<NASGU>
rng(4201,'twister');
nScenarios=5;
xiDraw=GDFM_DemoResidualBootstrap(R.idiosyncratic,24,nScenarios);
scenarios=zeros([size(X),nScenarios]);
for s=1:nScenarios
    Zscenario=R.common+xiDraw(:,:,s);
    scenarios(:,:,s)=bsxfun(@plus,bsxfun(@times,Zscenario,R.scale),R.mean);
end
stats=GDFM_DemoScenarioDiagnostics(X,scenarios);
R.exampleSensitivity=struct('parameter','Dynamic eigenvector count q', ...
    'qAlternative',qAlternative,'commonAlternative',Ralt.commonWithMean, ...
    'filterAlternative',Kalternative,'relativeCommonChange',commonDifference);
R.exampleScenarios=scenarios;
R.exampleScenarioStats=stats;
fprintf('  [TwoSided extension] changing q %d -> %d changes common by %.3g (relative Frobenius).\n', ...
    q,qAlternative,commonDifference);
fprintf('  Conditional block-bootstrap: %d scenarios; mean covariance error %.3f, lag-1 error %.3f.\n', ...
    nScenarios,mean(stats.covarianceRelativeError),mean(stats.meanAbsLag1CorrelationError));
if showPlots
    ix=1:min(250,size(X,1));
    figure('Name','TwoSided extension: filter sensitivity and conditional scenarios','Color','w');
    subplot(2,2,1);
    plot(ix,R.commonWithMean(ix,1),'b','LineWidth',1.1); hold on;
    plot(ix,Ralt.commonWithMean(ix,1),'r--','LineWidth',1.1);
    legend(sprintf('Original q=%d',q),sprintf('Alternative q=%d',qAlternative),'Location','best');
    xlabel('Time'); title('Changing retained dynamic subspace'); grid on;
    subplot(2,2,2);
    plot(R.filters.lags,klag,'b','LineWidth',1.1); hold on;
    plot(Ralt.filters.lags,altLagNorm,'r--','LineWidth',1.1);
    legend('Original K(L)','Alternative K(L)','Location','best');
    xlabel('Lag'); ylabel('Frobenius norm'); title('Dynamic projector filter sensitivity'); grid on;
    subplot(2,2,3);
    plot(ix,X(ix,1),'k-','LineWidth',1.1); hold on;
    plot(ix,R.commonWithMean(ix,1),'b--','LineWidth',1.1);
    plot(ix,scenarios(ix,1,1),'r');
    legend('Observed','Fixed common','Bootstrap scenario','Location','best');
    xlabel('Time'); title('Conditional scenario: series 1'); grid on;
    subplot(2,2,4);
    plot(ix,squeeze(scenarios(ix,1,:)),'LineWidth',0.8); hold on;
    plot(ix,R.commonWithMean(ix,1),'k','LineWidth',1.4);
    xlabel('Time'); title('Five scenarios around the same common path'); grid on;
end

end
