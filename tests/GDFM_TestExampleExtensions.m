function result = GDFM_TestExampleExtensions()
%GDFM_TESTEXAMPLEEXTENSIONS Native MATLAB smoke/identity tests for examples.
%
% Run from toolbox root after adding core/, examples/, tests/ to path:
%   addpath(genpath(pwd)); GDFM_TestExampleExtensions();
%
% Each example is executed without rendering figures. The four GDFM
% estimation algorithms remain untouched. These tests check dimensions,
% controlled changes, finite values, VAR recursion, PSD validity, and the
% Monte-Carlo covariance target. They do NOT establish scenario realism.

root=fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(fullfile(root,'examples'));
result=struct();

R=GDFM_Example_TwoSided(false,false);
assert(isfield(R,'exampleSensitivity') && R.exampleSensitivity.qAlternative~=R.settings.q);
assert(R.exampleSensitivity.relativeCommonChange>1e-8);
assert(size(R.exampleScenarios,3)==5 && all(isfinite(R.exampleScenarios(:))));
assert(norm(R.exampleScenarios(:,:,1)-R.exampleScenarios(:,:,2),'fro')>1e-6);
result.TwoSided=R.exampleSensitivity.relativeCommonChange;
fprintf('PASS: TwoSided projector sensitivity + conditional bootstrap\n');

R=GDFM_Example_OneSided(false,false);
assert(R.exampleSensitivity.relativeCommonChange>1e-8);
assert(all(isfinite(R.exampleScenarios(:))));
assert(size(R.exampleScenarioFactors,2)==R.settings.r);
assert(R.var.spectralRadius<1);
result.OneSided=R.exampleSensitivity.relativeCommonChange;
fprintf('PASS: OneSided factor intervention + finite-VAR generation\n');

R=GDFM_Example_Lee(false,false);
assert(R.exampleSensitivity.factorRecursionError<1e-8);
assert(norm(R.exampleSensitivity.commonImpulseEffect,'fro')>1e-6);
assert(all(isfinite(R.exampleScenarios(:))));
assert(size(R.exampleScenarioFactors,2)==R.settings.q);
result.Lee=R.exampleSensitivity.factorRecursionError;
fprintf('PASS: Lee shock propagation + Gaussian-shock generation\n');

R=GDFM_Example_Frequency(false,false);
S=R.exampleSensitivity.modifiedSpectrum;
assert(isempty(R.common)); % spectral decomposition is NOT time-domain common extraction
assert(norm(S(:)-R.spectral.full(:))>1e-6);
assert(R.exampleSensitivity.relativeFullRankProjectionError<1e-9);
assert(all(isfinite(R.exampleScenarios(:))));
for j=1:size(S,3)
    E=eig((S(:,:,j)+S(:,:,j)')/2);
    assert(min(real(E))>-1e-8);
end
assert(R.exampleSpectralSynthesis.relativeEnsembleCovarianceError<0.35);
result.Frequency=R.exampleSpectralSynthesis.relativeEnsembleCovarianceError;
fprintf('PASS: Frequency valid CPSD + distinct SPCA + spectral Gaussian synthesis\n');

fprintf('ALL GDFM example extension tests passed.\n');
end
