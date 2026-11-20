function R = GDFM_OneSidedFactors(X, spectral, r)
%GDFM_ONESIDEDFACTORS Weighted static-factor representation of common data.
%   Based on the spectral common covariance and generalized eigenvectors
%   of (Gamma_common(0), diagonal(Gamma_idiosyncratic(0))).
%
%   F = X*P and Lambda = Gamma_common(0)*P / (P'*Gamma_common(0)*P),
%   so common = F*Lambda' and idio = X-common.
%
%   This routine estimates contemporaneous in-sample static factors and
%   loadings; it does not implement the full forecasting algorithm of
%   Forni et al. [2]. Factors and loadings depend on the normalization
%   of the generalized eigenvectors.

[~,n] = size(X);
if r < 1 || r >= n || r ~= floor(r)
    error('GDFM:OneSided:r','r must be between 1 and n-1.');
end
C0 = real(mean(spectral.common,3));
C0 = (C0+C0')/2;
G0 = (X'*X)/size(X,1);
Xi0 = (G0-C0 + (G0-C0)')/2;

% Protect the generalized eigenproblem from zero/negative sample weights.
weightsRaw = real(diag(Xi0));
floorW = 1e-8*max(1,max(abs(diag(G0))));
weights = max(weightsRaw,floorW);
Winvhalf = diag(1./sqrt(weights));
H = Winvhalf*C0*Winvhalf;
H = (H+H')/2;
[U,D] = eig(H);
[evals,idx] = sort(real(diag(D)),'descend');
P = Winvhalf*U(:,idx(1:r));
Gram = P'*C0*P;
Gram = (Gram+Gram')/2;
% A pseudo-inverse is used if some selected directions have near-zero
% common power; a warning and the Gram matrix expose such cases.
if rcond(Gram) < 1e-10
    warning('GDFM:OneSided:Rank','Selected common directions are ill conditioned.');
end
Lambda = C0*P*pinv(Gram);
F = X*P;
common = F*Lambda';

R = struct();
R.common = real(common);
R.idio = X-R.common;
R.factors = real(F);
R.loadings = real(Lambda);
R.factorWeights = real(P);
R.genEigenvalues = evals(1:r);
R.commonCovariance0 = C0;
R.idioCovariance0 = Xi0;
R.weightsBeforeFloor = weightsRaw;
R.weightFloor = floorW;
R.normalizationGram = Gram;
end
