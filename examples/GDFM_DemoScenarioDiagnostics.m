function stats = GDFM_DemoScenarioDiagnostics(reference, scenarios)
%GDFM_DEMOSCENARIODIAGNOSTICS Descriptive, not pass/fail, scenario measures.
%   reference T-by-n; scenarios T-by-n-by-S in the SAME original units.
%   Compares empirical cross-series covariance and mean lag-1 correlation.
%   Neither statistic alone establishes adequate scenario realism.
[T,n] = size(reference);
if ~isnumeric(reference) || ~isreal(reference) || any(~isfinite(reference(:))) || ...
        ndims(scenarios)>3 || size(scenarios,1)~=T || ...
        size(scenarios,2)~=n || any(~isfinite(scenarios(:)))
    error('GDFM:DemoDiagnostics:Input','Reference/scenario dimensions or values are invalid.');
end
S = size(scenarios,3);
Cr = cov(reference,1);
lag1Reference = local_lag1(reference);
covError = zeros(S,1);
lag1Error = zeros(S,1);
for s = 1:S
    Ys = scenarios(:,:,s);
    covError(s) = norm(cov(Ys,1)-Cr,'fro')/max(norm(Cr,'fro'),eps);
    lag1Error(s) = mean(abs(local_lag1(Ys)-lag1Reference));
end
stats = struct('covarianceRelativeError',covError, ...
    'meanAbsLag1CorrelationError',lag1Error, ...
    'referenceLag1Correlation',lag1Reference);
end

function result = local_lag1(X)
[T,n] = size(X);
result = zeros(1,n);
for j = 1:n
    a = X(1:T-1,j); b = X(2:T,j);
    a = a-mean(a); b = b-mean(b);
    result(j) = (a'*b)/max(sqrt((a'*a)*(b'*b)),eps);
end
end
