function results = GDFM_CheckInstallation()
%GDFM_CHECKINSTALLATION Run the packaged MATLAB installation checks.
%
%   results = GDFM_CheckInstallation()
%
% Adds the release folders to the path and runs the core extraction,
% frequency-domain/SPCA, and scenario-example regression suites.
% The example tests use synthetic data and suppress plots.
% Run GDFM_RunAllExamples() separately to display demonstrations.

root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
fprintf('\nGDFM MATLAB Toolbox - installation checks\n');
fprintf('========================================\n');
GDFM_TestVariants();
GDFM_TestFrequencySPCA();
results = GDFM_TestExampleExtensions();
fprintf('\nInstallation checks completed successfully.\n');
end
