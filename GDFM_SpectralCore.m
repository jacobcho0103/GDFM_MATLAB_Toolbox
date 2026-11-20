function S = GDFM_SpectralCore(X, q, m, f, covmethod)
%GDFM_SPECTRALCORE Estimate the cross-spectrum and dynamic eigenstructure.
%   S = GDFM_SpectralCore(X,q,m,f,'Bartlett_triang')
%
%   X: T-by-n centered or standardized observations (stationarity assumed).
%   q: number of dynamic eigenvectors, 1 <= q < n.
%   m: covariance lag half-width, 1 <= m < T.
%   f: frequency-grid half-width; 2*f+1 equally spaced frequencies.
%
%   Returned fields (F = 2*f+1):
%     full, common, idio       n-by-n-by-F cross-spectral density matrices
%     eigvec                   n-by-q-by-F dynamic eigenvectors
%     eigval                   q-by-F eigenvalues in descending order
%     omega                    1-by-F angular frequency grid
%     covarianceLag            n-by-n-by-(2*m+1) weighted covariance sequence
%
%   Spectral factorization is S_common(w) = V_q(w)*D_q(w)*V_q(w)^H.
%   The frequency-domain *factor* in [4, eq. (9)], is
%   V_q(w)^H*S_full(w), which is a q-by-n matrix, NOT a time-series path.
%
%   A Bartlett covariance window is the recommended default. The 'Normal'
%   unweighted truncated covariance estimator can have negative eigenvalues.
%   They are retained rather than silently modified, to preserve the exact
%   identities A*B*S_full = S_common and S_full=S_common+S_idio.
%
%   References: [1] and [4] in README.md.

if nargin < 5, covmethod = 'Bartlett_triang'; end
[T,n] = size(X);
if ~isnumeric(X) || ~isreal(X) || ~ismatrix(X) || any(~isfinite(X(:)))
    error('GDFM:Spectral:Data','X must be a finite real matrix.');
end
if ~isscalar(q) || q ~= fix(q) || q < 1 || q >= n
    error('GDFM:Spectral:q','Require integer 1 <= q < n.');
end
if ~isscalar(m) || m ~= fix(m) || m < 1 || m >= T || ...
        ~isscalar(f) || f ~= fix(f) || f < 1
    error('GDFM:Spectral:Grid','Require 1 <= m < T and integer f >= 1.');
end
if ~(strcmpi(covmethod,'Bartlett_triang') || strcmpi(covmethod,'Normal'))
    error('GDFM:Spectral:Window','Use Bartlett_triang or Normal.');
end

Gamma = zeros(n,n,2*m+1);
for k = 0:m
    Gk = (X(k+1:T,:)' * X(1:T-k,:)) / (T-k);
    if strcmpi(covmethod,'Bartlett_triang')
        Gk = Gk*(m+1-k)/(m+1);
    end
    Gamma(:,:,m+1+k) = Gk;
    if k > 0, Gamma(:,:,m+1-k) = Gk'; end
end

nf = 2*f+1;
omega = 2*pi*(-f:f)/nf;
SX = complex(zeros(n,n,nf));
Sc = complex(zeros(n,n,nf));
V = complex(zeros(n,q,nf));
D = zeros(q,nf);

% Compute negative frequencies and omega=0, mirroring conjugate partners.
for j = 1:f+1
    Sj = complex(zeros(n,n));
    for k = -m:m
        Sj = Sj + Gamma(:,:,m+1+k)*exp(-1i*k*omega(j));
    end
    Sj = (Sj + Sj')/2;
    if j == f+1, Sj = real(Sj); end
    [U,Lambda] = eig(Sj);
    [sorted,order] = sort(real(diag(Lambda)),'descend');
    U = U(:,order(1:q));
    if j == f+1, U = real(U); end
    Sjcommon = U*diag(sorted(1:q))*U';
    SX(:,:,j) = Sj;
    V(:,:,j) = U;
    D(:,j) = sorted(1:q);
    Sc(:,:,j) = Sjcommon;
    if j <= f
        jm = nf+1-j;
        SX(:,:,jm) = conj(Sj);
        V(:,:,jm) = conj(U);
        D(:,jm) = D(:,j);
        Sc(:,:,jm) = conj(Sjcommon);
    end
end
S = struct('full',SX,'common',Sc,'idio',SX-Sc,'eigvec',V, ...
    'eigval',D,'omega',omega,'covarianceLag',Gamma, ...
    'lags',-m:m,'window',covmethod);
end
