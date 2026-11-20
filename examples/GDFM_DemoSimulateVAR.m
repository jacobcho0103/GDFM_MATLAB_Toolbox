function factors = GDFM_DemoSimulateVAR(V, innovations, initialFactors)
%GDFM_DEMOSIMULATEVAR Propagate a fitted finite VAR with user-provided noise.
%
% factors = GDFM_DemoSimulateVAR(V, innovations, initialFactors)
%   V               Result struct from GDFM_FactorVAR (coefficient,
%                   intercept, varLags, spectralRadius).
%   innovations     T-by-r finite matrix. First p rows are ignored.
%   initialFactors  p-by-r initial factor states in the SAME coordinates
%                   as V; here p=V.varLags. Can use V's fitted factor data.
%   factors         T-by-r simulated factor trajectories.
%
% f_t = c + sum_{j=1}^p M_j*f_(t-j) + epsilon_t.
% This is an additional finite-VAR model, not the GDFM estimator itself.
% A stationary process requires VAR stability and suitable innovations;
% drawing noise does not guarantee a realistic physical scenario.

if ~isstruct(V) || ~all(isfield(V,{'coefficient','intercept','varLags','spectralRadius'}))
    error('GDFM:DemoVAR:Model','V must be a fitted GDFM_FactorVAR result.');
end
p = V.varLags;
[T,r] = size(innovations);
if T <= p || ~isnumeric(innovations) || ~isreal(innovations) || ...
        any(~isfinite(innovations(:))) || ...
        size(V.coefficient,1)~=r || size(V.coefficient,2)~=r || ...
        size(V.coefficient,3)~=p || ...
        numel(V.intercept)~=r || ...
        ~isequal(size(initialFactors),[p,r]) || ...
        any(~isfinite(initialFactors(:)))
    error('GDFM:DemoVAR:Dimensions','Inconsistent factor, noise, or VAR dimensions.');
end
if V.spectralRadius >= 1
    error('GDFM:DemoVAR:Unstable','Cannot demonstrate stable VAR simulation (radius>=1).');
end
factors = zeros(T,r);
factors(1:p,:) = initialFactors;
for t = p+1:T
    value = V.intercept(:) + innovations(t,:)';
    for j = 1:p
        value = value + V.coefficient(:,:,j)*factors(t-j,:)';
    end
    factors(t,:) = value';
end
end
