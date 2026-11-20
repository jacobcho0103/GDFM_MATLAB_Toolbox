function result = GDFM_TestSeasonalPreprocessing()
%GDFM_TESTSEASONALPREPROCESSING Check optional MATLAB STL integration.
%
%   result = GDFM_TestSeasonalPreprocessing()
%
% Requires MATLAB R2021b or later for trenddecomp. Uses synthetic input,
% checks seasonal component dimensions and exact removal/restoration,
% and runs the standalone GDFM preprocessing demonstration without plots.

root=fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root,'examples'));
if verLessThan('matlab','9.11')
    error('GDFM:TestSeasonality:Version', ...
        'This optional test requires MATLAB R2021b or newer.');
end

T=720;
t=(0:T-1)';
x=[1.2*sin(2*pi*t/24)+0.3*cos(2*pi*t/60)+0.002*t, ...
   0.8*cos(2*pi*t/24)+0.5*sin(2*pi*t/60)-0.001*t, ...
   0.6*sin(2*pi*t/24)+0.4*cos(2*pi*t/60)+0.0005*t];

% Single-period STL: seasonality-only adjustment is the default.
[y1,p1]=GDFM_PreprocessSeasonality(x,24);
assert(~p1.removeTrend);
assert(isequal(size(p1.seasonal),size(x)));
assert(norm(GDFM_RestoreSeasonality(y1,p1)-x,'fro')<1e-10);

% Multiple-period STL with trend removal and multiple scenario paths.
[y2,p2]=GDFM_PreprocessSeasonality(x,[24 60],'RemoveTrend',true);
assert(p2.removeTrend && size(p2.seasonalByPeriod,3)==2);
assert(norm(p2.seasonal-sum(p2.seasonalByPeriod,3),'fro')<1e-10);
assert(norm(y2+p2.removed-x,'fro')<1e-10);
assert(norm(GDFM_RestoreSeasonality(y2,p2)-x,'fro')<1e-10);
scenarioSample=cat(3,y2,y2+0.1);
restored=GDFM_RestoreSeasonality(scenarioSample,p2);
assert(isequal(size(restored),[T 3 2]));
assert(norm(restored(:,:,1)-x,'fro')<1e-10);
assert(norm(restored(:,:,2)-(x+0.1),'fro')<1e-10);

% Standalone, no-figure GDFM example with an adjusted input.
E=GDFM_Example_Seasonality(false);
assert(E.restorationError<1e-10);
assert(isequal(size(E.scenariosRestored),[960 12 5]));
assert(all(isfinite(E.scenariosRestored(:))));
assert(norm(E.scenariosRestored(:,:,1)-E.scenariosRestored(:,:,2),'fro')>1e-6);
result=struct('singlePeriodError',norm(GDFM_RestoreSeasonality(y1,p1)-x,'fro'), ...
    'multiPeriodError',norm(GDFM_RestoreSeasonality(y2,p2)-x,'fro'), ...
    'exampleRestorationError',E.restorationError);
fprintf('PASS: STL seasonal adjustment, trend option, multiple periods, and scenario restoration\n');
end
