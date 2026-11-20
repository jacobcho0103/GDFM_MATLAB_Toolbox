function results = GDFM_RunAllExamples(varargin)
%GDFM_RUNALLEXAMPLES Run the GDFM examples and optional STL tutorial.
%
%   results = GDFM_RunAllExamples();
%   results = GDFM_RunAllExamples('IncludeSeasonality',true);
%
%   The optional STL demonstration needs MATLAB R2021b or later. The four
%   method-specific examples use stationary synthetic input and run
%   unchanged when IncludeSeasonality is false (default).
%
%   Folder structure:
%       GDFM_MATLAB_Toolbox/
%           GDFM_Extract.m, GDFM_SpectralCore.m, ... (core)
%           examples/
%               GDFM_DemoData.m
%               GDFM_SelectQ.m
%               GDFM_Example_TwoSided.m
%               GDFM_Example_OneSided.m
%               GDFM_Example_Lee.m
%               GDFM_Example_Frequency.m
%               GDFM_Example_Seasonality.m (optional, R2021b+)
%               GDFM_RunAllExamples.m
%
%   No additional MATLAB toolboxes or external datasets are required.
%   GDFM_PreprocessSeasonality is a separate, optional preprocessing step.
%   The four method examples cover factor selection, decomposition,
%   sensitivity, and illustrative scenarios. The STL example shows
%   seasonal adjustment and restoration around a GDFM fit.

p=inputParser;
p.FunctionName='GDFM_RunAllExamples';
addParameter(p,'IncludeSeasonality',false, ...
    @(v) (islogical(v)||isnumeric(v)) && isscalar(v) && (v==0 || v==1));
parse(p,varargin{:});

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
if p.Results.IncludeSeasonality
    results.Seasonality=GDFM_Example_Seasonality();
end
fprintf('\nGDFM demonstrations finished. Inspect the figures.\n');
fprintf('Scenario diagnostics are included in each returned result.\n');
end
