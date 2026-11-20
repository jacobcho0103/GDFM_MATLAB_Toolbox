function R = GDFM_FactorVAR(F,p,identifyShocks,estimator)
%GDFM_FACTORVAR Estimate the factor VAR and optionally Cholesky shocks.
%   V = GDFM_FactorVAR(F,p,true,'YuleWalker')
%   V = GDFM_FactorVAR(F,p,false,'OLS')
%
%   Factor model (column-vector convention):
%      f_t = c + M_1 f_(t-1) + ... + M_p f_(t-p) + epsilon_t,
%      C(L) = I - M_1 L - ... - M_p L^p.
%   With Cholesky identification: epsilon_t=Z delta_t, E[delta delta']=I.
%
%   'YuleWalker' implements covariance block equations (14)-(16) in
%   Lee & Baldick (2017), IEEE TPWRS 32(1), for a finite VAR(p).
%   'OLS' retains the previous least-squares estimator for compatibility.
%   Both use the fitted factor mean as a constant/intercept; the factors
%   themselves are NOT changed. First p innovation rows are NaN.
%
%   NOTE: The white-noise/shock interpretation is an estimation ASSUMPTION.
%   Cholesky alone whitens the contemporaneous covariance and does not
%   establish time independence or unique structural shock identification.
%   Factor rotation, ordering, and normalization affect M, Z, and delta.

if nargin<3, identifyShocks=false; end
if nargin<4 || isempty(estimator), estimator='OLS'; end
[T,r] = size(F);
if ~isnumeric(F) || ~isreal(F) || isempty(F) || any(~isfinite(F(:)))
    error('GDFM:VAR:Data','F must be finite and real.');
end
if ~isscalar(p) || ~isfinite(p) || p~=fix(p) || p<1 || T-p<=r*p
    error('GDFM:VAR:Lags','Require integer p>=1 and T-p > r*p.');
end
if ~(strcmpi(estimator,'OLS') || strcmpi(estimator,'YuleWalker'))
    error('GDFM:VAR:Estimator','Estimator must be OLS or YuleWalker.');
end
mu = mean(F,1);
Fc = bsxfun(@minus,F,mu);
Y = Fc(p+1:T,:);
regressors = zeros(T-p,r*p);
for j = 1:p
    regressors(:,(j-1)*r+(1:r)) = Fc(p+1-j:T-j,:);
end
if rank(regressors) < r*p
    error('GDFM:VAR:Rank','Lagged factors are rank deficient.');
end
momentError = NaN;
if strcmpi(estimator,'OLS')
    coef = regressors\Y;
    M = zeros(r,r,p);
    for j = 1:p
        M(:,:,j) = coef((j-1)*r+(1:r),:)';
    end
else
    % Biased lag-covariance sequence Gamma(h)=(1/T)sum_t f_t f_(t-h)^T.
    % Constant denominator yields a positive-semidefinite block covariance.
    Gamma = zeros(r,r,p+1);
    for h = 0:p
        Gamma(:,:,h+1) = (Fc(h+1:T,:)'*Fc(1:T-h,:))/T;
    end
    H = zeros(r,r*p);
    G = zeros(r*p,r*p);
    for h = 1:p
        H(:,(h-1)*r+(1:r)) = Gamma(:,:,h+1);
        for j = 1:p
            d = h-j;
            if d>=0
                block = Gamma(:,:,d+1);
            else
                block = Gamma(:,:,-d+1)';
            end
            G((j-1)*r+(1:r),(h-1)*r+(1:r)) = block;
        end
    end
    if rcond(G) < 1e-12
        error('GDFM:VAR:YWConditioning', ...
          'Yule-Walker covariance block matrix is nearly singular.');
    end
    Mstack = H/G;
    momentError = norm(Mstack*G-H,'fro')/max(1,norm(H,'fro'));
    M = zeros(r,r,p);
    for j = 1:p
        M(:,:,j) = Mstack(:,(j-1)*r+(1:r));
    end
end
C = zeros(r,r,p+1);
C(:,:,1) = eye(r);
for j = 1:p, C(:,:,j+1) = -M(:,:,j); end

% Recover intercept c=(I-sum_j M_j)*mean(f) as a row vector.
Msum = zeros(r);
for j = 1:p, Msum = Msum + M(:,:,j); end
intercept = (eye(r)-Msum)*mu';
fit = NaN(T,r);
fit(p+1:T,:) = repmat(intercept',T-p,1);
for j = 1:p
    fit(p+1:T,:) = fit(p+1:T,:) + F(p+1-j:T-j,:)*M(:,:,j)';
end
eps = NaN(T,r);
eps(p+1:T,:) = F(p+1:T,:) - fit(p+1:T,:);
valid = eps(p+1:T,:);
innovCov = (valid'*valid)/(T-p);
innovCov = (innovCov+innovCov')/2;

companion = zeros(r*p);
companion(1:r,:) = reshape(M,r,r*p);
if p>1, companion(r+1:r*p,1:r*(p-1))=eye(r*(p-1)); end
spectralRadius = max(abs(eig(companion)));
if spectralRadius>=1
    warning('GDFM:VAR:Unstable', ...
      'Estimated VAR spectral radius %.3g >= 1.',spectralRadius);
end

shockLoading = [];
whiteNoise = [];
if identifyShocks
    [shockLoading,flag] = chol(innovCov,'lower');
    if flag~=0 || rcond(shockLoading)<1e-12
        error('GDFM:VAR:Covariance', ...
          'Innovations covariance is not positive definite.');
    end
    whiteNoise=NaN(T,r);
    whiteNoise(p+1:T,:)=(shockLoading\valid')';
end
R = struct('coefficient',M,'noiseLoading',C,'intercept',intercept', ...
    'innovations',eps,'shockLoading',shockLoading, ...
    'whitenedShocks',whiteNoise,'innovationCovariance',innovCov, ...
    'fitted',fit,'spectralRadius',spectralRadius,'varLags',p, ...
    'innovationStart',p+1,'estimator',estimator, ...
    'yuleWalkerMomentError',momentError);
end
