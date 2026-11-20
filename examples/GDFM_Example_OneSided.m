function R = GDFM_Example_OneSided(showScreening,showPlots)
%GDFM_EXAMPLE_ONESIDED Example 2: contemporaneous one-sided factor estimate.
%
%   R = GDFM_Example_OneSided();
%   R = GDFM_Example_OneSided(false); % omit repeated q-screening figure
%   R = GDFM_Example_OneSided(false,false); % compute without figures
%
%   Demonstrates q-screening, static factor extraction and loadings,
%   and common/idiosyncratic decomposition, with an optional finite OLS VAR.
%   IMPORTANT: Current OneSided branch implements the contemporaneous
%   in-sample estimator associated with Forni (2005); it is NOT the
%   complete one-sided causal forecasting algorithm.
%   q is the dynamic eigenvector count; r is the number of static factors.
%   Here we choose r=q for ease of interpretation, not as a general rule.

%   Scenario example: intervene on a static factor, then generate scenarios by
%   simulating the additional fitted OLS VAR with resampled innovations.

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
m=16; f=40;
[q,screen] = GDFM_SelectQ(X,'QMax',5,'m',m,'f',f,'Plot',showScreening); %#ok<NASGU>
r = q;
R = GDFM_Extract(X,'Method','OneSided','q',q,'r',r, ...
    'm',m,'f',f,'VARLags',2);

trueCommon = bsxfun(@plus,truth.common,truth.mean);
relativeError = norm(R.commonWithMean-trueCommon,'fro')/norm(truth.common,'fro');
fprintf('\n[OneSided] q=%d, r=%d (synthetic truth q=%d)\n',q,r,truth.q);
fprintf('  Relative error vs known common: %.4f\n',relativeError);
fprintf('  Static factors: %d x %d; static loadings: %d x %d\n', ...
    size(R.factors),size(R.loadings));
fprintf('  f_t = (X_standardized)_t * P; chi_t = f_t * Lambda^T.\n');
fprintf('  Generalized eigenvalues (selected): '); fprintf('%.3g ',R.static.genEigenvalues); fprintf('\n');
fprintf('  VAR(p) is an additional modeling step, not required for decomposition.\n');

ix=1:min(250,size(X,1));
if showPlots
figure('Name','Example 2: OneSided - static factors and loadings','Color','w');
subplot(2,3,1);
plot(ix,X(ix,1),'Color',[0.65 0.65 0.65]); hold on;
plot(ix,trueCommon(ix,1),'k--','LineWidth',1.2);
plot(ix,R.commonWithMean(ix,1),'b','LineWidth',1.1);
legend('Observed','True common','Estimated common','Location','best');
xlabel('Time index'); title('One-sided common estimate'); grid on;
subplot(2,3,2);
plot(ix,truth.idiosyncratic(ix,1),'k--'); hold on;
plot(ix,R.idiosyncraticOriginal(ix,1),'r');
legend('True residual','Estimated residual','Location','best');
xlabel('Time index'); title('Idiosyncratic component'); grid on;
subplot(2,3,3);
plot(ix,R.factors(ix,1:min(2,r)));
xlabel('Time index'); title(sprintf('Static factors f_t (%d extracted)',r)); grid on;
subplot(2,3,4);
imagesc(R.loadings); colorbar;
xlabel('Factor index'); ylabel('Observed series'); title('Static loadings Lambda');
subplot(2,3,5);
bar(1:r,R.static.genEigenvalues);
xlabel('Static factor index'); ylabel('Generalized eigenvalue');
title('Generalized eigenstructure'); grid on;
subplot(2,3,6);
valid=R.var.innovationStart:size(X,1);
plot(valid,R.innovations(valid,1),'Color',[0.25 0.45 0.8]);
xlabel('Time index'); ylabel('VAR residual');
title('Optional VAR innovations: first factor'); grid on;
end % decomposition figures

%% Scenario and sensitivity demonstration: factor intervention and finite-VAR scenario simulation
% Controlled factor perturbation: keep static A0 and the observed residual
% fixed. This demonstrates factor-loadings transmission, not probability.
Fchanged=R.factors;
Fchanged(:,1)=1.30*Fchanged(:,1);
changedCommon=Fchanged*R.loadings';
changedOriginal=bsxfun(@plus,bsxfun(@times, ...
    changedCommon+R.idiosyncratic,R.scale),R.mean);
interventionChange=norm(changedCommon-R.common,'fro')/max(1,norm(R.common,'fro'));

% Additional scenario model: fitted finite OLS VAR on static factor scores.
% Joint-row innovation resampling retains within-time factor correlation;
% idiosyncratic residuals use moving-block resampling. Both assume a
% fixed fitted model; neither establishes full Forni causal forecasting.
if R.var.spectralRadius>=1
    error('GDFM:Example:UnstableVAR', ...
        'OneSided demo VAR unstable; choose a different fitted VAR before simulation.');
end
oldState=rng;
cleanup=onCleanup(@() rng(oldState)); %#ok<NASGU>
rng(4202,'twister');
[T,n]=size(X);
p=R.var.varLags;
epsPool=R.innovations(p+1:end,:);
nScenarios=5;
xiDraw=GDFM_DemoResidualBootstrap(R.idiosyncratic,24,nScenarios);
scenarios=zeros(T,n,nScenarios);
factorScenarios=zeros(T,r,nScenarios);
for s=1:nScenarios
    epsDraw=zeros(T,r);
    epsDraw(p+1:T,:)=epsPool(randi(size(epsPool,1),T-p,1),:);
    factorScenarios(:,:,s)=GDFM_DemoSimulateVAR(R.var,epsDraw,R.factors(1:p,:));
    Zscenario=factorScenarios(:,:,s)*R.loadings'+xiDraw(:,:,s);
    scenarios(:,:,s)=bsxfun(@plus,bsxfun(@times,Zscenario,R.scale),R.mean);
end
stats=GDFM_DemoScenarioDiagnostics(X,scenarios);
R.exampleSensitivity=struct('factorGain',1.30, ...
    'factorIntervention',Fchanged,'intervenedObservations',changedOriginal, ...
    'relativeCommonChange',interventionChange);
R.exampleScenarios=scenarios;
R.exampleScenarioFactors=factorScenarios;
R.exampleScenarioStats=stats;
fprintf('  [OneSided extension] factor 1 multiplied by 1.30; common change %.3g.\n',interventionChange);
fprintf('  Finite OLS-VAR + residual blocks: %d scenarios; mean covariance error %.3f, lag-1 error %.3f.\n', ...
    nScenarios,mean(stats.covarianceRelativeError),mean(stats.meanAbsLag1CorrelationError));
if showPlots
    ix=1:min(250,T);
    figure('Name','OneSided extension: factor intervention and finite-VAR scenarios','Color','w');
    subplot(2,2,1);
    plot(ix,R.factors(ix,1),'b'); hold on;
    plot(ix,Fchanged(ix,1),'r--');
    legend('Original factor 1','Factor 1 x 1.30','Location','best');
    xlabel('Time'); title('Controlled static-factor intervention'); grid on;
    subplot(2,2,2);
    plot(ix,X(ix,1),'k'); hold on;
    plot(ix,changedOriginal(ix,1),'r');
    legend('Observed','Changed factor + same residual','Location','best');
    xlabel('Time'); title('Propagation through static loadings A_0'); grid on;
    subplot(2,2,3);
    plot(ix,R.factors(ix,1),'k','LineWidth',1.1); hold on;
    plot(ix,squeeze(factorScenarios(ix,1,:)),'LineWidth',0.8);
    xlabel('Time'); title('Simulated factor trajectories from VAR'); grid on;
    subplot(2,2,4);
    plot(ix,X(ix,1),'k','LineWidth',1.3); hold on;
    plot(ix,squeeze(scenarios(ix,1,:)),'LineWidth',0.8);
    xlabel('Time'); title('Finite-VAR conditional model scenarios'); grid on;
end

end
