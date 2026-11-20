function [q, info] = GDFM_SelectQ(X, varargin)
%GDFM_SELECTQ Exploratory dynamic-eigenvalue gap screening for q.
%
%   [q,info] = GDFM_SelectQ(X,'QMax',5,'m',16,'f',40,'Plot',true)
%
%   Estimates the CPSD of standardized X, averages its ordered dynamic
%   eigenvalues over frequencies, and chooses the largest adjacent ratio
%             lambda_j / lambda_(j+1),    j = 1,...,QMax.
%   This is a transparent SCREE/GAP HEURISTIC, not an optimality theorem,
%   a K-means factor-count selection procedure, or a formal information
%   criterion. Check the plot and use domain knowledge before fixing q.
%   For OneSided, the static factor count r can differ from q; examples
%   set r=q for clarity, not because it is mathematically required.
%
%   Inputs: X must be T-by-n, with at least four nonconstant columns.
%   Output info.eigenvalues (n-1-by-1), info.ratios (QMax-by-1),
%   info.candidates, info.settings, and info.methodDescription.
%
%   Uses GDFM_SpectralCore.m (base MATLAB; no Statistics Toolbox required).

[T,n] = size(X);
if ~isnumeric(X) || ~isreal(X) || ~ismatrix(X) || any(~isfinite(X(:))) || n<4 || T<20
    error('GDFM:SelectQ:Data','Provide a finite real T-by-n matrix (T>=20, n>=4).');
end
p = inputParser;
addParameter(p,'QMax',min(5,n-2),@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v==fix(v)&&v>=1&&v<=n-2);
addParameter(p,'m',min(16,floor(sqrt(T))),@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v==fix(v)&&v>=1&&v<T);
addParameter(p,'f',40,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v==fix(v)&&v>=1);
addParameter(p,'Plot',true,@(v) islogical(v)&&isscalar(v));
parse(p,varargin{:});
op = p.Results;

sx = std(X,0,1);
if any(sx<=0)
    error('GDFM:SelectQ:Constant','Remove constant columns first.');
end
Xz = bsxfun(@rdivide,bsxfun(@minus,X,mean(X,1)),sx);
S = GDFM_SpectralCore(Xz,n-1,op.m,op.f,'Bartlett_triang');
% Clamp only for the DIAGNOSTIC, to avoid divide-by-zero from roundoff.
lam = mean(max(real(S.eigval),0),2);
base = max(1,max(lam));
lamFloor = max(1e-12*base,eps(base));
candidates = (1:op.QMax)';
ratios = lam(candidates)./max(lam(candidates+1),lamFloor);
[~,idx] = max(ratios);
q = candidates(idx);
info = struct('eigenvalues',lam,'ratios',ratios, ...
    'candidates',candidates,'recommendedQ',q, ...
    'settings',struct('m',op.m,'f',op.f,'QMax',op.QMax), ...
    'methodDescription','Largest adjacent mean dynamic-eigenvalue ratio; heuristic only');

fprintf('GDFM factor-number screening: suggested q = %d (largest mean-eigenvalue ratio %.3g).\n', ...
    q,ratios(idx));
fprintf('Inspect the spectral-gap plot when choosing q for your data.\n');
if op.Plot
    figure('Name','GDFM factor-number screening','Color','w');
    subplot(1,2,1);
    semilogy(1:numel(lam),max(lam,lamFloor),'o-','LineWidth',1.2);
    hold on; plot(q, max(lam(q),lamFloor),'rs','MarkerSize',8,'LineWidth',1.5);
    xlabel('Dynamic eigenvalue index'); ylabel('Mean eigenvalue (log scale)');
    title('Mean CPSD dynamic eigenspectrum'); grid on;
    subplot(1,2,2);
    plot(candidates,ratios,'o-','LineWidth',1.2); hold on;
    plot(q,ratios(idx),'rs','MarkerSize',8,'LineWidth',1.5);
    xlabel('Candidate factor count q'); ylabel('lambda_q / lambda_{q+1}');
    title(sprintf('Spectral gap: suggested q = %d',q)); grid on;
end
end
