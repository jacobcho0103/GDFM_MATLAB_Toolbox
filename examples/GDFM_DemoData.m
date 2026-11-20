function [X, truth] = GDFM_DemoData(varargin)
%GDFM_DEMODATA Reproducible stationary multivariate time series with known factors.
%
%   [X,truth] = GDFM_DemoData();
%   [X,truth] = GDFM_DemoData('T',960,'n',12,'Seed',2026,'NoiseStd',0.25);
%
%   Creates n observations driven by two dynamic factors with AR(2)
%   dynamics, contemporaneous/lagged loadings and idiosyncratic AR(1)
%   noise. This is simulated statistical data, NOT wind power or a
%   physically constrained power-production process.
%
%   X: T-by-n observed time series (including nonzero means).
%   truth.common: T-by-n common component, excluding the mean.
%   truth.idiosyncratic: T-by-n idiosyncratic component.
%   truth.factors: T-by-2 latent factor trajectories.
%   truth.loading0/loading1: n-by-2 loading matrices, lag 0 and 1.
%   truth.factorInnovations: T-by-2 factor driving noises.
%   truth.mean: 1-by-n deterministic offsets.
%
%   Since factor coordinates are rotation/sign ambiguous, do not expect
%   the estimated factor columns to equal truth.factors column-by-column.
%   The generator restores the caller's random-number state.

p = inputParser;
p.FunctionName = 'GDFM_DemoData';
addParameter(p,'T',960,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v==fix(v)&&v>=100);
addParameter(p,'n',12,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v==fix(v)&&v>=6);
addParameter(p,'Seed',2026,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v==fix(v)&&v>=0);
addParameter(p,'NoiseStd',0.25,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0);
parse(p,varargin{:});
T = p.Results.T;
n = p.Results.n;
oldState = rng;
cleanup = onCleanup(@() rng(oldState)); %#ok<NASGU>
rng(p.Results.Seed,'twister');

% Two common AR(2) factors with distinct temporal behaviors.
u = randn(T,2);
F = zeros(T,2);
for t = 3:T
    F(t,1) = 0.76*F(t-1,1) - 0.10*F(t-2,1) + u(t,1);
    F(t,2) = -0.42*F(t-1,2) + 0.18*F(t-2,2) + u(t,2);
end

% Distinct spatial loading patterns; the lagged term creates dynamics.
i = (1:n)';
A0 = [0.95 + 0.24*sin(2*pi*i/n), 0.95*cos(2*pi*i/n)];
A1 = [0.22 + 0.08*cos(2*pi*i/n), 0.18*sin(2*pi*i/n)];
chi = F*A0';
chi(2:end,:) = chi(2:end,:) + F(1:end-1,:)*A1';

% Local AR(1) idiosyncratic processes (uncorrelated across locations).
xi = zeros(T,n);
noise = p.Results.NoiseStd*randn(T,n);
for t = 2:T
    xi(t,:) = 0.25*xi(t-1,:) + noise(t,:);
end

baseline = 1 + 0.3*cos(2*pi*i'/n);
X = chi + xi + repmat(baseline,T,1);
truth = struct('common',chi,'idiosyncratic',xi,'factors',F, ...
    'factorInnovations',u,'loading0',A0,'loading1',A1, ...
    'mean',baseline,'q',2,'T',T,'n',n, ...
    'Seed',p.Results.Seed,'NoiseStd',p.Results.NoiseStd);
end
