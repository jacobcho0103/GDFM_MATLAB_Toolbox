function [Xadjusted, prep] = GDFM_PreprocessSeasonality(X, periods, varargin)
%GDFM_PREPROCESSSEASONALITY Remove seasonal structure using MATLAB STL.
%
%   [Xadjusted,prep] = GDFM_PreprocessSeasonality(X,24)
%   [Xadjusted,prep] = GDFM_PreprocessSeasonality(X,[24 168], ...
%       'RemoveTrend',true)
%
%   Applies MATLAB's trenddecomp(...,'stl',periods) separately to each
%   column of a T-by-n matrix. Periods are positive integer numbers of
%   equally spaced observations (for hourly data, 24 means one day).
%   STL is provided by MATLAB R2021b and newer; no separate toolbox is
%   required. This is OPTIONAL preprocessing, not part of GDFM estimation.
%
%   Default: remove the seasonal estimate, but retain the estimated trend.
%   'RemoveTrend',true also removes the estimated long-term trend, which
%   can be useful when preparing approximately stationary GDFM inputs.
%
%   prep.trend            : T-by-n estimated long-term trend
%   prep.seasonal         : T-by-n sum of seasonal estimates
%   prep.seasonalByPeriod : T-by-n-by-numel(periods) seasonal estimates
%   prep.remainder        : T-by-n STL remainder
%   prep.removed          : T-by-n deterministic part subtracted
%   prep.periods          : STL periods, in observations
%   prep.removeTrend      : whether the trend is also removed
%
%   Recover the input with GDFM_RestoreSeasonality(Xadjusted,prep).
%   That utility also restores same-horizon T-by-n-by-S scenario arrays.
%   Estimation/restoration of future seasonal paths is outside the scope
%   of this in-sample utility; future periods require extrapolation.
%
%   See also trenddecomp, GDFM_RestoreSeasonality, GDFM_Example_Seasonality.

if ~(isnumeric(X) && ismatrix(X) && isreal(X) && ~isempty(X) && ...
        all(isfinite(X(:))))
    error('GDFM:Seasonality:Input', ...
        'X must be a finite, real, nonempty T-by-n numeric matrix.');
end
if ~(isnumeric(periods) && isvector(periods) && ~isempty(periods) && ...
        isreal(periods) && all(isfinite(periods(:))) && ...
        all(periods(:)==fix(periods(:))) && all(periods(:)>=2) && ...
        numel(unique(periods(:)))==numel(periods))
    error('GDFM:Seasonality:Periods', ...
        'periods must be distinct integers of at least 2 observations.');
end
periods=double(periods(:)');
[T,n]=size(X);
if T < 3*max(periods)
    error('GDFM:Seasonality:Length', ...
        'Supply at least three complete cycles of the longest STL period.');
end
p=inputParser;
p.FunctionName='GDFM_PreprocessSeasonality';
addParameter(p,'RemoveTrend',false, ...
    @(v) (islogical(v) || isnumeric(v)) && isscalar(v) && ...
          isfinite(double(v)) && (v==0 || v==1));
parse(p,varargin{:});
removeTrend=logical(p.Results.RemoveTrend);

if verLessThan('matlab','9.11')
    error('GDFM:Seasonality:Version', ...
        'STL preprocessing requires MATLAB R2021b or newer (trenddecomp).');
end

X=double(X);
nPeriods=numel(periods);
longterm=zeros(T,n);
seasonalByPeriod=zeros(T,n,nPeriods);
remainder=zeros(T,n);
for j=1:n
    [lt,st,r]=trenddecomp(X(:,j),'stl',periods);
    st=reshape(st,T,[]);
    if size(st,2)~=nPeriods
        error('GDFM:Seasonality:STLOutput', ...
            'STL returned an unexpected number of seasonal components.');
    end
    longterm(:,j)=lt(:);
    seasonalByPeriod(:,j,:)=reshape(st,T,1,nPeriods);
    remainder(:,j)=r(:);
end
seasonal=sum(seasonalByPeriod,3);
removed=seasonal;
if removeTrend
    removed=removed+longterm;
end
Xadjusted=X-removed;
if ~all(isfinite(Xadjusted(:)))
    error('GDFM:Seasonality:Nonfinite','STL generated nonfinite values.');
end
prep=struct('method','STL','periods',periods, ...
    'removeTrend',removeTrend,'trend',longterm, ...
    'seasonal',seasonal,'seasonalByPeriod',seasonalByPeriod, ...
    'remainder',remainder,'removed',removed, ...
    'inputSize',[T n]);
end
