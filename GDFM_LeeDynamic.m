function R = GDFM_LeeDynamic(X, spectral, q, lags)
%GDFM_LEEDYNAMIC Extract dynamic factors and loadings from spectral DPCA.
%   R = GDFM_LeeDynamic(X,S,q,-m:m)
%
%   Implements the finite-filter extraction in [3, eqs. (2), (3), (10), (13)].
%   Given spectral common-space eigenvectors V_q(omega):
%      A(omega) = V_q(omega), B(omega) = V_q(omega)^H
%      f_t = B(L) X_t, chi_t = A(L) f_t.
%   Inverse frequency transforms obtain finite-lag A(L), B(L).
%   Positive lags are past observations, negative lags are leads.
%
%   INPUTS
%     X        T-by-n centered/standardized observations.
%     spectral Structure from GDFM_SpectralCore with eigvec and omega.
%     q        Number of dynamic factors (columns in V_q).
%     lags     Integer vector of retained filter lags, e.g. -m:m.
%
%   OUTPUTS
%     factors         T-by-q time-domain dynamic factor scores B(L)X.
%     common          T-by-n estimated common component A(L)f.
%     loadings        n-by-q-by-numel(lags) finite impulse response A(L).
%     extraction      q-by-n-by-numel(lags) finite impulse response B(L).
%     spectralLoading n-by-q-by-F phase-aligned spectral eigenvectors.
%     spectralFilter  q-by-n-by-F spectral extraction filter.
%     spectrumCommon  n-by-n-by-F spectral common component.
%     phaseInvarianceError  projection difference before/after alignment.
%
%   IDENTIFICATION/APPROXIMATION:
%   Dynamic eigenvectors are defined only up to frequency-specific phase
%   and (for repeated eigenvalues) unitary rotation. Here adjacent positive
%   frequencies are aligned column-wise by phase, with conjugate symmetry
%   imposed for negative frequencies. This changes the extracted factor
%   coordinates, but not the common spectral projection. Near eigenvalue
%   crossings, more general rotation ambiguity remains.
%   Truncating A and B separately and zero-padding boundaries means that
%   A(L)B(L)X generally differs from the direct two-sided projector K(L)X.
%   The identity chi=A(L)f is exact for the RETURNED finite representation,
%   not proof of asymptotic GDFM estimation consistency.
%
%   No scenario generation or forecasting is performed.

[T,n] = size(X);
[nSpectral,qSpectral,nf] = size(spectral.eigvec);
if nSpectral ~= n || q ~= qSpectral || numel(spectral.omega) ~= nf
    error('GDFM:Lee:Dimensions','Spectral dimensions differ from X or q.');
end
if any(lags ~= fix(lags)) || numel(unique(lags)) ~= numel(lags)
    error('GDFM:Lee:Lags','Lags must be distinct integers.');
end
if nf ~= 2*floor(nf/2)+1
    error('GDFM:Lee:Grid','The frequency grid must have odd length.');
end
f = (nf-1)/2;
A = spectral.eigvec;
% Start from the real eigenvectors at zero frequency. Phase-align at
% increasing positive frequencies, then mirror by conjugation.
for j = f+2:nf
    for k = 1:q
        overlap = A(:,k,j-1)'*A(:,k,j);
        if abs(overlap) > 1e-12
            A(:,k,j) = A(:,k,j)*exp(-1i*angle(overlap));
        end
    end
    A(:,:,nf+1-j) = conj(A(:,:,j));
end
B = complex(zeros(q,n,nf));
Scommon = complex(zeros(n,n,nf));
maxProjDiff = 0;
for j = 1:nf
    B(:,:,j) = A(:,:,j)';
    projector = A(:,:,j)*B(:,:,j);
    baseline = spectral.eigvec(:,:,j)*spectral.eigvec(:,:,j)';
    maxProjDiff = max(maxProjDiff,norm(projector-baseline,'fro'));
    Scommon(:,:,j) = projector*spectral.full(:,:,j);
end
H_A = local_ifourier(A,spectral.omega,lags);
H_B = local_ifourier(B,spectral.omega,lags);
F = real(GDFM_ApplyFilter(X,H_B,lags));
chi = real(GDFM_ApplyFilter(F,H_A,lags));

R = struct();
R.factors = F;
R.common = chi;
R.idio = X-chi;
R.loadings = H_A;
R.extraction = H_B;
R.spectralLoading = A;
R.spectralFilter = B;
R.spectrumCommon = Scommon;
R.lags = lags;
R.phaseInvarianceError = maxProjDiff;
R.timeSamples = T;
end

function H = local_ifourier(Hspec,omega,lags)
[a,b,nf] = size(Hspec);
H = zeros(a,b,numel(lags));
for k = 1:numel(lags)
    tmp = complex(zeros(a,b));
    for j = 1:nf
        tmp = tmp + Hspec(:,:,j)*exp(1i*lags(k)*omega(j));
    end
    H(:,:,k) = real(tmp/nf);
end
end
