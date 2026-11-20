function samples = GDFM_DemoResidualBootstrap(residual, blockLength, nScenarios)
%GDFM_DEMORESIDUALBOOTSTRAP Joint circular moving-block residual resampling.
%
% samples = GDFM_DemoResidualBootstrap(residual,blockLength,nScenarios)
%   residual     T-by-n real finite estimated idiosyncratic component.
%   blockLength  Positive integer <= T. Larger blocks preserve longer local
%                temporal dependence. Does not preserve dynamics exactly.
%   nScenarios   Number of independent block-bootstrap scenarios.
%   samples      T-by-n-by-nScenarios resampled paths in input units.
%
% Entire ROW VECTORS are resampled together, preserving within-block
% cross-series structure. Residual columns are centered before sampling.
% The supplied residual pool and its block boundaries are not a complete
% stochastic model. Bootstrapped residuals can exhibit discontinuities at
% block boundaries; they need not be independent of the common component.
% Caller controls rng for reproducibility. Base MATLAB only.

if ~isnumeric(residual) || ~isreal(residual) || ~ismatrix(residual) || ...
        isempty(residual) || any(~isfinite(residual(:)))
    error('GDFM:DemoBootstrap:Input','Residual input must be finite real T-by-n.');
end
[T,n] = size(residual);
if ~isscalar(blockLength) || blockLength ~= fix(blockLength) || ...
        blockLength < 1 || blockLength > T
    error('GDFM:DemoBootstrap:Block','blockLength must be integer in [1,T].');
end
if ~isscalar(nScenarios) || nScenarios ~= fix(nScenarios) || nScenarios < 1
    error('GDFM:DemoBootstrap:Count','nScenarios must be positive integer.');
end
pool = bsxfun(@minus,residual,mean(residual,1));
samples = zeros(T,n,nScenarios);
for s = 1:nScenarios
    indices = zeros(T,1);
    pos = 1;
    while pos <= T
        start = randi(T);
        count = min(blockLength,T-pos+1);
        indices(pos:pos+count-1) = mod(start-1+(0:count-1),T)+1;
        pos = pos+count;
    end
    samples(:,:,s) = pool(indices,:);
end
end
