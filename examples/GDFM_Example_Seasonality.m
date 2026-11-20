function result = GDFM_Example_Seasonality(showPlots)
%GDFM_EXAMPLE_SEASONALITY Demonstrate optional STL before GDFM fitting.
%
%   result = GDFM_Example_Seasonality();      % with figures
%   result = GDFM_Example_Seasonality(false); % without figures
%
%   Adds a known 24-observation seasonal profile and a small deterministic
%   trend to the stationary GDFM synthetic data. Uses MATLAB trenddecomp
%   through GDFM_PreprocessSeasonality to remove them before fitting Lee's
%   GDFM formulation. Demonstrates restoration after residual sampling.
%   Any of the four GDFM methods can be applied to the adjusted data.
%
%   Requires MATLAB R2021b or later for the optional STL utility.

if nargin<1, showPlots=true; end
if ~(isscalar(showPlots) && (islogical(showPlots) || ...
        (isnumeric(showPlots) && (showPlots==0 || showPlots==1))))
    error('GDFM:Example:PlotFlag','showPlots must be a logical scalar.');
end
exampleDir=fileparts(mfilename('fullpath'));
addpath(exampleDir);
addpath(fileparts(exampleDir));

[Xbase,truth]=GDFM_DemoData('Seed',2026);
[T,n]=size(Xbase);
t=(0:T-1)';
spatial=(1:n);
% Known additive periodic pattern: the phase and magnitude vary by series.
seasonalTrue=sin(2*pi*t/24)*(0.8+0.2*cos(2*pi*spatial/n)) + ...
    cos(2*pi*t/24)*(0.3*sin(2*pi*spatial/n));
trendTrue=((t-(T-1)/2)/T)*(0.4+0.1*cos(2*pi*spatial/n));
Xraw=Xbase+seasonalTrue+trendTrue;

% STL is an independent preprocessing step, like factor-number screening.
[Xadjusted,prep]=GDFM_PreprocessSeasonality(Xraw,24,'RemoveTrend',true);
R=GDFM_Extract(Xadjusted,'Method','Lee','q',truth.q,'r',truth.q, ...
    'm',16,'f',40,'VARLags',2);

% Example conditional scenarios: keep Lee's fitted common component,
% bootstrap joint residual blocks, convert back to adjusted-data units,
% then restore the estimated deterministic seasonal/trend components.
oldState=rng;
cleanup=onCleanup(@() rng(oldState)); %#ok<NASGU>
rng(4301,'twister');
nScenarios=5;
xiDraw=GDFM_DemoResidualBootstrap(R.idiosyncratic,24,nScenarios);
scenariosAdjusted=zeros(T,n,nScenarios);
for s=1:nScenarios
    standardized=R.common+xiDraw(:,:,s);
    scenariosAdjusted(:,:,s)=bsxfun(@plus, ...
        bsxfun(@times,standardized,R.scale),R.mean);
end
scenariosRestored=GDFM_RestoreSeasonality(scenariosAdjusted,prep);
recovered=GDFM_RestoreSeasonality(Xadjusted,prep);
reconstructionError=norm(recovered-Xraw,'fro')/max(1,norm(Xraw,'fro'));

result=struct('original',Xraw,'base',Xbase,'knownSeasonal',seasonalTrue, ...
    'knownTrend',trendTrue,'adjusted',Xadjusted, ...
    'preprocessing',prep,'model',R, ...
    'scenariosAdjusted',scenariosAdjusted, ...
    'scenariosRestored',scenariosRestored, ...
    'restorationError',reconstructionError);
fprintf('\n[Seasonality] MATLAB STL: period=24, trend removal enabled.\n');
fprintf('  Input and adjusted observations: %d time steps, %d series.\n',T,n);
fprintf('  In-sample seasonal/trend restoration error: %.3e.\n',reconstructionError);
fprintf('  Lee GDFM fit on adjusted data; %d scenarios restored to observed scale.\n',nScenarios);
if showPlots
    idx=1:min(240,T);
    figure('Name','GDFM preprocessing: MATLAB STL seasonal adjustment','Color','w');
    subplot(2,2,1);
    plot(idx,Xraw(idx,1),'Color',[0.5 0.5 0.5]); hold on;
    plot(idx,Xadjusted(idx,1),'b');
    legend('Raw data','Adjusted data','Location','best');
    title('Observed vs seasonally adjusted series 1');
    xlabel('Time index'); grid on;
    subplot(2,2,2);
    plot(idx,seasonalTrue(idx,1),'k--'); hold on;
    plot(idx,prep.seasonal(idx,1),'b');
    legend('Known seasonal','STL estimate','Location','best');
    title('Seasonal component'); xlabel('Time index'); grid on;
    subplot(2,2,3);
    plot(t,trendTrue(:,1)+truth.mean(1),'k--'); hold on;
    plot(t,prep.trend(:,1),'r');
    legend('Known baseline + trend','STL estimate','Location','best');
    title('Long-term trend (series 1)'); xlabel('Time index'); grid on;
    subplot(2,2,4);
    plot(idx,Xraw(idx,1),'k','LineWidth',1.2); hold on;
    plot(idx,squeeze(scenariosRestored(idx,1,:)),'LineWidth',0.8);
    title('Scenarios after restoring estimated seasonality');
    xlabel('Time index'); grid on;
end
end
