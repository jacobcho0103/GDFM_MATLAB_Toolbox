function results = GDFM_CheckInstallation()
%GDFM_CHECKINSTALLATION Run the packaged MATLAB installation checks.
%
%   results = GDFM_CheckInstallation()
%
% Adds the release folders to the path and runs the core extraction,
% frequency-domain/SPCA, and scenario-example regression suites.
% Optional STL seasonal-preprocessing tests run on MATLAB R2021b or newer.
% The example tests use synthetic data and suppress plots.
% Run GDFM_RunAllExamples() separately to display demonstrations.

root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
fprintf('\nGDFM MATLAB Toolbox - installation checks\n');
fprintf('========================================\n');
GDFM_TestVariants();
GDFM_TestFrequencySPCA();
results = GDFM_TestExampleExtensions();
if verLessThan('matlab','9.11')
    fprintf('SKIP: STL seasonal preprocessing (requires MATLAB R2021b or newer).\n');
else
    results.Seasonality = GDFM_TestSeasonalPreprocessing();
end
fprintf('\nInstallation checks completed successfully.\n');
end
