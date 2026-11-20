function results = GDFM_RunAllExamples()
%GDFM_RUNALLEXAMPLES Execute all four standalone GDFM demonstration files.
%
%   results = GDFM_RunAllExamples();
%
%   Folder structure:
%       GDFM_MATLAB_Toolbox_v2.0.0/
%           GDFM_Extract.m, GDFM_SpectralCore.m, ... (core)
%           examples/
%               GDFM_DemoData.m
%               GDFM_SelectQ.m
%               GDFM_Example_TwoSided.m
%               GDFM_Example_OneSided.m
%               GDFM_Example_Lee.m
%               GDFM_Example_Frequency.m
%               GDFM_RunAllExamples.m
%
%   No additional MATLAB toolboxes or external datasets are required.
%   Each example covers factor selection, decomposition, controlled
%   sensitivity analysis, and an illustrative scenario-generation model.

examplesDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(examplesDir);
addpath(rootDir);
addpath(examplesDir);
assert(exist('GDFM_Extract','file')==2, ...
    'GDFM core source files are missing from the parent folder.');

fprintf('\n=== GDFM decomposition, sensitivity, and scenario examples ===\n');
results = struct();
results.TwoSided  = GDFM_Example_TwoSided();
results.OneSided  = GDFM_Example_OneSided(false);
results.Lee       = GDFM_Example_Lee(false);
results.Frequency = GDFM_Example_Frequency(false);
fprintf('\nAll four demonstrations finished. Inspect the original and extension figures.\n');
fprintf('Scenario diagnostics are included in each returned result.\n');
end
