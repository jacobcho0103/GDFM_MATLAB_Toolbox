function Xrestored = GDFM_RestoreSeasonality(Xadjusted, prep)
%GDFM_RESTORESEASONALITY Restore in-sample STL seasonal/trend components.
%
%   Xrestored = GDFM_RestoreSeasonality(Xadjusted,prep)
%
%   Xadjusted can be T-by-n or T-by-n-by-S (S scenarios). The T-by-n
%   deterministic component prep.removed, obtained from
%   GDFM_PreprocessSeasonality, is added to every scenario.
%
%   Only the sample horizon used to fit STL is supported. To generate a
%   future horizon, construct an explicit future seasonal/trend forecast
%   rather than reusing the in-sample estimated components.
%
%   See also GDFM_PreprocessSeasonality.

if ~isstruct(prep) || ~isfield(prep,'removed') || ...
        ~isfield(prep,'inputSize')
    error('GDFM:Seasonality:Info', ...
        'prep must be returned by GDFM_PreprocessSeasonality.');
end
if ~(isnumeric(Xadjusted) && isreal(Xadjusted) && ...
        all(isfinite(Xadjusted(:))) && ndims(Xadjusted)<=3)
    error('GDFM:Seasonality:Input', ...
        'Xadjusted must be a finite, real T-by-n or T-by-n-by-S array.');
end
[T,n]=size(prep.removed);
if size(Xadjusted,1)~=T || size(Xadjusted,2)~=n || ...
        ~isequal(prep.inputSize,[T n])
    error('GDFM:Seasonality:Size', ...
        'Time and series dimensions must match the fitted STL sample.');
end
Xrestored=bsxfun(@plus,double(Xadjusted),reshape(prep.removed,T,n,1));
end
