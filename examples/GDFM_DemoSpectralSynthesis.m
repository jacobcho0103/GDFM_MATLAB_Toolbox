function [scenarios, info] = GDFM_DemoSpectralSynthesis(S, omega, T, nScenarios)
%GDFM_DEMOSPECTRALSYNTHESIS Illustrative stationary Gaussian CPSD sampler.
%
% [Y,info] = GDFM_DemoSpectralSynthesis(S,omega,T,nScenarios)
%   S          n-by-n-by-F Hermitian PSD CPSD, symmetric in frequencies.
%   omega      Equally spaced angular frequency grid -pi..pi, centered
%              on zero as used by GDFM_SpectralCore (F odd).
%   T          Even positive integer >= 4, desired time-series length.
%   nScenarios Number of independent Gaussian realizations.
%   Y          T-by-n-by-nScenarios, zero-mean stationary Gaussian paths.
%   info       Target zero-lag covariance after spectrum interpolation,
%              and mean empirical covariance across scenarios.
%
% Procedure: interpolate the complex Hermitian CPSD onto FFT frequencies,
% project each interpolated matrix onto the PSD cone (clipping tiny
% negative eigenvalues), generate conjugate-symmetric complex Gaussian
% DFT coefficients with E[dft*dft^H]=T*S(omega), and apply IFFT.
% This generates zero-mean Gaussian samples whose EXPECTED finite-grid
% second-order structure approximates the supplied CPSD. A single finite
% realization will not match the spectrum exactly. Interpolation, periodic
% DFT boundary assumptions and finite samples introduce discrepancies.
% Caller controls rng. This is a downstream DEMONSTRATION, not GDFM fitting,
% the GAN-based method of Cho et al., or a validated wind-power model.

[n,n2,F] = size(S);
if ~isnumeric(S) || n~=n2 || mod(F,2)~=1 || F<3 || ...
        any(~isfinite(S(:))) || ~isscalar(T) || T~=fix(T) || ...
        mod(T,2)~=0 || T<4 || ~isscalar(nScenarios) || ...
        nScenarios~=fix(nScenarios) || nScenarios<1
    error('GDFM:DemoSpectral:Input','Invalid spectrum or even T/scenario count.');
end
f=(F-1)/2;
wExpected=2*pi*(-f:f)/F;
if numel(omega)~=F || norm(omega(:)'-wExpected,inf)>1e-9
    error('GDFM:DemoSpectral:Grid','omega must be 2*pi*(-f:f)/(2*f+1).');
end
scale=max(1,norm(S(:)));
for j=1:F
    if norm(S(:,:,j)-S(:,:,j)','fro')>1e-8*scale || ...
            norm(S(:,:,j)-conj(S(:,:,F+1-j)),'fro')>1e-8*scale
        error('GDFM:DemoSpectral:Symmetry','Spectrum is not Hermitian/conjugate symmetric.');
    end
end
% Precompute one positive-frequency square root at each FFT bin.
nbin=T/2+1;
wTarget=2*pi*(0:T/2)/T;
wPos=omega(f+1:end);
rootSpec=complex(zeros(n,n,nbin));
targetCov=zeros(n);
for k=1:nbin
    w=min(wTarget(k),wPos(end));
    hi=find(wPos>=w,1,'first');
    if isempty(hi), hi=numel(wPos); end
    if hi==1
        M=S(:,:,f+1);
    else
        lo=hi-1;
        alpha=(w-wPos(lo))/(wPos(hi)-wPos(lo));
        M=(1-alpha)*S(:,:,f+lo)+alpha*S(:,:,f+hi);
    end
    M=(M+M')/2;
    if k==1 || k==nbin
        M=real(M);
        M=(M+M')/2;
    end
    [V,D]=eig(M);
    d=real(diag(D));
    if min(d)<-1e-7*max(1,norm(M,2))
        error('GDFM:DemoSpectral:PSD','Target spectrum has a negative eigenvalue.');
    end
    d=max(d,0);
    rootSpec(:,:,k)=V*diag(sqrt(d));
    Mcorrected=V*diag(d)*V';
    weight=2;
    if k==1 || k==nbin, weight=1; end
    targetCov=targetCov+weight*real(Mcorrected)/T;
end
scenarios=zeros(T,n,nScenarios);
meanCov=zeros(n);
for s=1:nScenarios
    Fourier=complex(zeros(T,n));
    for k=1:nbin
        if k==1 || k==nbin
            z=randn(n,1);
        else
            z=(randn(n,1)+1i*randn(n,1))/sqrt(2);
        end
        y=sqrt(T)*rootSpec(:,:,k)*z;
        Fourier(k,:)=y.';
        if k>1 && k<nbin
            Fourier(T-k+2,:)=conj(y.');
        end
    end
    path=real(ifft(Fourier,[],1));
    scenarios(:,:,s)=path;
    meanCov=meanCov+(path'*path)/T/nScenarios;
end
info=struct('targetCovariance',targetCov, ...
    'meanEmpiricalCovariance',meanCov, ...
    'relativeEnsembleCovarianceError',norm(meanCov-targetCov,'fro')/max(norm(targetCov,'fro'),eps), ...
    'frequencyBins',nbin);
end
