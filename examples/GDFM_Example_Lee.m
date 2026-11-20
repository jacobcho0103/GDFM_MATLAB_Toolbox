function R = GDFM_Example_Lee(showScreening,showPlots)
%GDFM_EXAMPLE_LEE Example 3: Lee-Baldick (2017) dynamic factor extraction.
%
%   R = GDFM_Example_Lee();
%   R = GDFM_Example_Lee(false); % omit repeated q-screening figure
%   R = GDFM_Example_Lee(false,false); % compute without figures
%
%   Uses q dynamic eigenvectors, lagged loading A(L), extraction B(L),
%   extracted factors f_t=B(L)X_t, and chi_t=A(L)f_t. Fits a finite
%   Yule-Walker VAR(p), extracts epsilon_t, and identifies shocks delta_t
%   using Cholesky epsilon_t=Z*delta_t.
%
%   NOTE: A and B are finite two-sided filters in this implementation,
%   not causal forecasts. Phase, truncation and Cholesky ordering affect
%   factor/shock coordinates. Initial p innovation/shock rows are NaN.

%   Scenario example: perturb an identified shock, propagate the factor VAR and
%   dynamic loading, then generate Gaussian-shock scenario paths.

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
m=16; f=40; p=2;
[q,screen] = GDFM_SelectQ(X,'QMax',5,'m',m,'f',f,'Plot',showScreening); %#ok<NASGU>
R = GDFM_Extract(X,'Method','Lee','q',q,'r',q, ...
    'm',m,'f',f,'VARLags',p);

trueCommon = bsxfun(@plus,truth.common,truth.mean);
relativeError = norm(R.commonWithMean-trueCommon,'fro')/norm(truth.common,'fro');
valid=R.var.innovationStart:size(X,1);
EPS=R.innovations(valid,:);
DEL=R.shocks(valid,:);
Erec=(R.shockLoading*DEL')';
shockIdentityError=norm(EPS-Erec,'fro')/max(1,norm(EPS,'fro'));
fprintf('\n[Lee] q=%d (synthetic truth q=%d)\n',q,truth.q);
fprintf('  Relative error vs known common: %.4f\n',relativeError);
fprintf('  Dynamic factors: %d x %d; A(L): %d x %d x %d\n', ...
    size(R.factors),size(R.loadings));
fprintf('  B(L): %d x %d x %d\n',size(R.filters.factorExtraction));
fprintf('  VAR estimation: %s; stability radius: %.3f\n', ...
    R.var.estimator,R.var.spectralRadius);
fprintf('  epsilon_t = Z*delta_t relative identity error: %.2e\n',shockIdentityError);
fprintf('  Delta_t is covariance-whitened, not proven temporally independent.\n');

ix=1:min(250,size(X,1));
if showPlots
figure('Name','Example 3: Lee - factors, filters, shocks','Color','w');
subplot(2,3,1);
plot(ix,X(ix,1),'Color',[0.65 0.65 0.65]); hold on;
plot(ix,trueCommon(ix,1),'k--','LineWidth',1.1);
plot(ix,R.commonWithMean(ix,1),'b','LineWidth',1.1);
legend('Observed','True common','Estimated common','Location','best');
xlabel('Time index'); title('Lee common estimate'); grid on;
subplot(2,3,2);
plot(ix,truth.idiosyncratic(ix,1),'k--'); hold on;
plot(ix,R.idiosyncraticOriginal(ix,1),'r');
legend('True residual','Estimated residual','Location','best');
xlabel('Time index'); title('Idiosyncratic component'); grid on;
subplot(2,3,3);
plot(ix,R.factors(ix,1:min(2,q)));
xlabel('Time index'); title('Dynamic factor trajectories B(L)X'); grid on;
subplot(2,3,4);
imagesc(R.filters.lags,1:q,reshape(R.loadings(1,:,:),[q,numel(R.filters.lags)]));
xlabel('Lag (negative = future)'); ylabel('Factor');
title('Dynamic loadings A(L): series 1'); colorbar;
subplot(2,3,5);
imagesc(R.filters.lags,1:q,reshape(R.filters.factorExtraction(:,1,:),[q,numel(R.filters.lags)]));
xlabel('Lag (negative = future)'); ylabel('Factor');
title('Dynamic filter B(L): series 1'); colorbar;
subplot(2,3,6);
plot(valid,R.innovations(valid,1),'k-'); hold on;
plot(valid,R.shocks(valid,1),'b-');
legend('Innovations epsilon','Whitened shocks delta','Location','best');
xlabel('Time index'); title('VAR innovation and identified shock'); grid on;
end % decomposition figures

%% Scenario and sensitivity demonstration: a shock intervention and shock-driven stochastic scenarios
% In the fitted VAR representation: epsilon=Z*delta, then
% f_t = intercept + sum_j M_j*f_(t-j) + epsilon_t and chi=A(L)*f.
% A one-time shock intervention propagates through the fitted dynamics.
[T,n]=size(X);
impulseTime=80;
shockGain=3.0;
baseEps=zeros(T,q);
baseEps(valid,:)=R.innovations(valid,:);
perturbedEps=baseEps;
perturbedEps(impulseTime,:)=perturbedEps(impulseTime,:)+ ...
    (R.shockLoading*[shockGain;zeros(q-1,1)])';
fBase=GDFM_DemoSimulateVAR(R.var,baseEps,R.factors(1:p,:));
fPert=GDFM_DemoSimulateVAR(R.var,perturbedEps,R.factors(1:p,:));
commonBase=real(GDFM_ApplyFilter(fBase,R.loadings,R.filters.lags));
commonPert=real(GDFM_ApplyFilter(fPert,R.loadings,R.filters.lags));
impulseCommonChange=commonPert-commonBase;

% Scenario-generation extension: draw new unit-covariance Gaussian shocks
% delta, map them through the estimated Cholesky Z, fit the factor VAR,
% apply A(L), and add jointly block-resampled idiosyncratic components.
% Gaussian shock independence and stationarity are model assumptions.
if R.var.spectralRadius>=1
    error('GDFM:Example:UnstableVAR', ...
        'Lee demo VAR unstable; cannot run stationary shock simulation.');
end
oldState=rng;
cleanup=onCleanup(@() rng(oldState)); %#ok<NASGU>
rng(4203,'twister');
nScenarios=5;
xiDraw=GDFM_DemoResidualBootstrap(R.idiosyncratic,24,nScenarios);
scenarios=zeros(T,n,nScenarios);
factorScenarios=zeros(T,q,nScenarios);
for s=1:nScenarios
    delta=zeros(T,q);
    delta(p+1:T,:)=randn(T-p,q);
    epsDraw=delta*R.shockLoading';
    factorScenarios(:,:,s)=GDFM_DemoSimulateVAR(R.var,epsDraw,R.factors(1:p,:));
    Zscenario=real(GDFM_ApplyFilter( ...
        factorScenarios(:,:,s),R.loadings,R.filters.lags))+xiDraw(:,:,s);
    scenarios(:,:,s)=bsxfun(@plus,bsxfun(@times,Zscenario,R.scale),R.mean);
end
stats=GDFM_DemoScenarioDiagnostics(X,scenarios);
R.exampleSensitivity=struct('shockIndex',1,'shockTime',impulseTime, ...
    'shockMagnitude',shockGain,'factorsPerturbed',fPert, ...
    'commonImpulseEffect',impulseCommonChange, ...
    'factorRecursionError',norm(fBase-R.factors,'fro')/max(1,norm(R.factors,'fro')));
R.exampleScenarios=scenarios;
R.exampleScenarioFactors=factorScenarios;
R.exampleScenarioStats=stats;
fprintf('  [Lee extension] shock 1 impulse %.1f at t=%d; common response norm %.3f.\n', ...
    shockGain,impulseTime,norm(impulseCommonChange,'fro'));
fprintf('  Fitted-VAR factor recursion relative error: %.3e.\n',R.exampleSensitivity.factorRecursionError);
fprintf('  Gaussian shocks + A(L) + residual blocks: %d scenarios; mean covariance error %.3f, lag-1 error %.3f.\n', ...
    nScenarios,mean(stats.covarianceRelativeError),mean(stats.meanAbsLag1CorrelationError));
if showPlots
    ix=1:min(250,T);
    figure('Name','Lee extension: shock propagation and scenario synthesis','Color','w');
    subplot(2,2,1);
    plot(ix,fBase(ix,1),'k'); hold on;
    plot(ix,fPert(ix,1),'r--','LineWidth',1.1);
    legend('Original extracted factor','With shock impulse','Location','best');
    xlabel('Time'); title('Propagation through finite factor VAR'); grid on;
    subplot(2,2,2);
    plot(ix,impulseCommonChange(ix,1),'b','LineWidth',1.2); hold on;
    plot([impulseTime impulseTime],ylim,'r--');
    xlabel('Time'); title('Effect on series 1 through dynamic A(L)'); grid on;
    subplot(2,2,3);
    plot(ix,R.factors(ix,1),'k','LineWidth',1.2); hold on;
    plot(ix,squeeze(factorScenarios(ix,1,:)),'LineWidth',0.7);
    xlabel('Time'); title('Generated factor 1 (black: extracted factor)'); grid on;
    subplot(2,2,4);
    plot(ix,X(ix,1),'k','LineWidth',1.2); hold on;
    plot(ix,squeeze(scenarios(ix,1,:)),'LineWidth',0.8);
    xlabel('Time'); title('Shock-driven illustrative scenarios'); grid on;
end

end
