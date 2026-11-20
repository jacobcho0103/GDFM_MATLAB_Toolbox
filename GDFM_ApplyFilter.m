function Y = GDFM_ApplyFilter(X, H, lags)
%GDFM_APPLYFILTER Apply a matrix-valued lag filter with zero boundary padding.
%   Y = GDFM_APPLYFILTER(X,H,LAGS) evaluates
%       Y(t,:)' = sum_k H(:,:,k) * X(t-LAGS(k),:)'
%   for X (T-by-nInput) and H (nOutput-by-nInput-by-nLag).
%   Positive LAGS use past observations; negative LAGS use future data.
%   Samples outside [1,T] are treated as zero. H can be complex, although
%   conjugate-symmetric spectral constructions should yield real H.

[T,n] = size(X);
[nOut,nIn,nLag] = size(H);
if nIn ~= n || numel(lags) ~= nLag
    error('GDFM:Filter:Dimensions','Filter/input dimensions are incompatible.');
end
Y = zeros(T,nOut);
for j = 1:nLag
    ell = lags(j);
    if ell ~= floor(ell)
        error('GDFM:Filter:Lags','All lags must be integers.');
    end
    t = max(1,1+ell):min(T,T+ell);
    if ~isempty(t)
        Y(t,:) = Y(t,:) + X(t-ell,:)*H(:,:,j).';
    end
end
end
