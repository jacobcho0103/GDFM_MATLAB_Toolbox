function O = GDFM_FrequencySPCA(SX, omega, Xstd, r)
%GDFM_FREQUENCYSPCA Zero-lag CPSD-to-covariance static-PCA projection.
%
% O = GDFM_FrequencySPCA(SX,omega,Xstd,r)
%
% Input
%   SX    n-by-n-by-F complex cross-power spectral density on the equally
%         spaced omega = 2*pi*(-f:f)/F grid, F = 2*f+1.
%   omega 1-by-F angular frequencies ordered from negative to positive.
%   Xstd  T-by-n centered (preferably standardized) observation data.
%   r     1 <= r <= n: number of retained STATIC eigenvectors. This is
%         independent of the DYNAMIC factor count q used to build SX.
%
% Procedure [4, eqs. (14)-(15), in dimensionally consistent form]:
%   1. IFFT SX to obtain covariance sequence Phi(k).
%   2. Select zero lag Phi(0), eigendecompose Phi(0) = P*Lambda*P'.
%   3. Obtain projected signal Xhat = (P_r*P_r')*X in n-by-T notation,
%      equivalently Xhat = Xstd*P_r*P_r' in T-by-n notation.
%      For positive selected eigenvalues, this equals
%      Phi(0)*P_r*inv(Lambda_r)*P_r'*X (n-by-T notation).
%
% IMPORTANT: This is a deterministic STATIC PROJECTION of observed data.
%   r=n recovers Xstd (even after modifying SX), because P*P'=I.
%   Changing eigenvalues alone does NOT change this projection.
%   Generating new realizations that match an altered CPSD requires a
%   separate synthesis procedure. A CPSD alone does not determine a
%   unique time-series realization.
%
% Outputs: covarianceLags (n-by-n-by-F), lags (-f:f), covariance0,
%   eigenvectors, eigenvalues, retainedEigenvectors, staticScores,
%   projectedStandardized, residualStandardized, projection,
%   relativeProjectionError, inverseFormValid.
%
% The code uses P_r*P_r' instead of an explicit matrix inverse for
% numerical stability and to handle rank-deficient spectra transparently.
%
% Requires MATLAB base functionality only. No scenario generation.

if ~isnumeric(SX) || ndims(SX) ~= 3 || any(~isfinite(SX(:)))
    error('GDFM:SPCA:Spectrum','SX must be a finite n-by-n-by-F array.');
end
[n,n2,F] = size(SX);
if n ~= n2 || mod(F,2) ~= 1 || F < 3
    error('GDFM:SPCA:Dimensions','SX must be square on an odd frequency grid.');
end
if ~isnumeric(Xstd) || ~isreal(Xstd) || size(Xstd,2) ~= n || ...
        any(~isfinite(Xstd(:)))
    error('GDFM:SPCA:Data','Xstd must be finite real T-by-n observations.');
end
if ~isscalar(r) || r ~= fix(r) || r < 1 || r > n
    error('GDFM:SPCA:Count','Static factor count must satisfy 1 <= r <= n.');
end
if ~isnumeric(omega) || numel(omega) ~= F || any(~isfinite(omega(:)))
    error('GDFM:SPCA:Grid','omega must contain F finite frequencies.');
end
f = (F-1)/2;
expected = 2*pi*(-f:f)/F;
if norm(omega(:).'-expected,inf) > 1e-9
    error('GDFM:SPCA:Grid','omega must equal 2*pi*(-f:f)/(2*f+1).');
end

% Check that SX is a real-signal CPSD (Hermitian and conjugate pairs).
tol = 1e-7 * max(1,norm(SX(:)));
for j = 1:F
    if norm(SX(:,:,j)-SX(:,:,j)','fro') > tol
        error('GDFM:SPCA:Hermitian','SX is not Hermitian at a frequency.');
    end
    if norm(SX(:,:,j)-conj(SX(:,:,F+1-j)),'fro') > tol
        error('GDFM:SPCA:Conjugate','SX lacks conjugate frequency symmetry.');
    end
end

% ifftshift puts the zero frequency in the first position for IFFT;
% fftshift returns covariance lags in order -f,...,f.
PhiComplex = fftshift(ifft(ifftshift(SX,3),[],3),3);
if max(abs(imag(PhiComplex(:)))) > tol
    error('GDFM:SPCA:Complex','IFFT covariance has non-negligible imaginary values.');
end
Phi = real(PhiComplex);
Phi0 = (Phi(:,:,f+1) + Phi(:,:,f+1)')/2;
[Q,D] = eig(Phi0);
[eigvals,order] = sort(real(diag(D)),'descend');
Q = real(Q(:,order));

% A negative eigenvalue means the proposed zero-lag covariance is invalid.
eigTol = 1e-9*max(1,norm(Phi0,2));
if min(eigvals) < -eigTol
    error('GDFM:SPCA:PSD','Zero-lag covariance is not positive semidefinite.');
end
P = Q(:,1:r);
Pj = P*P';
Xhat = Xstd*Pj;
resid = Xstd-Xhat;
O = struct();
O.covarianceLags = Phi;
O.lags = -f:f;
O.covariance0 = Phi0;
O.eigenvectors = Q;
O.eigenvalues = eigvals;
O.retainedEigenvectors = P;
O.staticScores = Xstd*P;
O.projection = Pj;
O.projectedStandardized = Xhat;
O.residualStandardized = resid;
O.relativeProjectionError = norm(resid,'fro')/max(norm(Xstd,'fro'),eps);
O.inverseFormValid = all(eigvals(1:r)>eigTol);
O.staticFactors = r;
end
