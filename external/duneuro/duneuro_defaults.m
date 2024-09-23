function cfg = duneuro_defaults(cfg)

% This function is based on the bst-duneuro functon bst_load_default_duneuroConfiguration
 
if nargin<1
  cfg = [];
end

%% General settings
cfg.outputpath                     = ft_getopt(cfg, 'outputpath', tempdir);
cfg.duneuro_configuration_filename = ft_getopt(cfg, 'duneuro_configuration_filename', 'duneuro_minifile.mini');
cfg.dnFemMethodType                = ft_getopt(cfg, 'dnFemMethodType',                'fitted'); % 'fitted' or 'unfitted'
cfg.dnFemSolverType                = ft_getopt(cfg, 'dnFemSolverType',                'cg');     % what else?
cfg.dnMeshElementType              = ft_getopt(cfg, 'dnMeshElementType',              []);       % this should be determined from the mesh
cfg.dnGeometryAdapted              = ft_getopt(cfg, 'dnGeometryAdapted',              false);  % true or  false why and  how
cfg.dnTolerance                    = ft_getopt(cfg, 'dnTolerance',                    1e-8);
cfg.EnableCacheMemory              = ft_getopt(cfg, 'EnableCacheMemory',              false); % no idea

%% 1 Sensors
% subpart electrode : [electrodes]
cfg.dnElectrodType                 = ft_getopt(cfg, 'dnElectrodType',                 'normal'); % FIXME: typo

%% subpart [meg]
cfg.dnMegIntorderadd               = ft_getopt(cfg, 'dnMegIntorderadd',               0); % FIXME: why is this numeric, is it intended boolean? why also false/true as strings?
cfg.dnMegType                      = ft_getopt(cfg, 'dnMegType',                      'physical');

%% 4 - Subpart  [solver] ==> refers to the linear system solver ?
cfg.dnSolverSolverType = ft_getopt(cfg, 'dnSolverSolverType',                         'cg'); %  conjugate solver  what are the others 
if ~isfield(cfg,'dnSolverPreconditionerType'); cfg.dnSolverPreconditionerType ='amg'; end %  what are the others 
if ~isfield(cfg,'dnSolverCgSmootherType'); cfg.dnSolverCgSmootherType ='ssor'; end %  what are the others 
if ~isfield(cfg,'dnSolverIntorderadd'); cfg.dnSolverIntorderadd =0; end %  what are the others 

% case of the dg discontinious galerkin
if ~isfield(cfg,'DgSmootherType'); cfg.dnGgSmootherType = 'ssor'; end %  what are the others 
if ~isfield(cfg,'DgScheme'); cfg.dnScheme = 'sipg'; end %  what are the others 
if ~isfield(cfg,'DgPenalty'); cfg.dnPenalty = 20; end %  what are the others 
if ~isfield(cfg,'DgEdgeNormType'); cfg.dnEdgeNormType = 'houston'; end %  what are the others 
if ~isfield(cfg,'DgWeights'); cfg.dnWeights = true; end %  what are the others 
if ~isfield(cfg,'DgReduction'); cfg.dnWeights = true; end %  what are the others 

%% 5 - Subpart  [solution]
if ~isfield(cfg,'dnSolutionPostProcess'); cfg.dnSolutionPostProcess = true; end %  what are the others 
if ~isfield(cfg,'dnSolutionSubstractMean'); cfg.dnSolutionSubstractMean =false; end %  what are the others 
% subpart  [solution.solver]
if ~isfield(cfg,'dnSolutionSolverReduction'); cfg.dnSolutionSolverReduction = 1e-10; end %  what are the others 

%% 6 - subpart  [solution.source_model]
if ~isfield(cfg,'femSourceModel'); cfg.femSourceModel = 'venant'; end % partial_integration, venant, subtraction | expand smtype
if ~isfield(cfg,'femSourceModelIntorderadd'); cfg.femSourceModelIntorderadd = 0; end % 
if ~isfield(cfg,'femSourceModelIntorderadd_lb'); cfg.femSourceModelIntorderadd_lb = 2; end % 
if ~isfield(cfg,'femSourceModelNumberOfMoments'); cfg.femSourceModelNumberOfMoments = 3; end % 
if ~isfield(cfg,'femSourceModelReferenceLength'); cfg.femSourceModelReferenceLength = 20; end % 
if ~isfield(cfg,'femSourceModelWeightingExponent'); cfg.femSourceModelWeightingExponent = 1; end % 
if ~isfield(cfg,'femSourceModelRelaxationFactor'); cfg.femSourceModelRelaxationFactor = 6; end % will be converted to 10^-(#)% 
if ~isfield(cfg,'femSourceModelMixedMoments'); cfg.femSourceModelMixedMoments = true; end % 
if ~isfield(cfg,'femSourceModelRestrict'); cfg.femSourceModelRestrict = true; end % 
if ~isfield(cfg,'femSourceModelInitialization'); cfg.femSourceModelInitialization = 'closest_vertex'; end % 

%% 7 - subpart  [brainstorm]
%if ~isfield(cfg,'brainstormModality'); cfg.brainstormModality = cfg.modality; end % should be done from the bst as tmp folder 
if ~isfield(cfg,'BstOutputFolder'); cfg.BstOutputFolder = cfg.outputpath; end % should be done from the bst as tmp folder 
if ~isfield(cfg,'BstSaveTransfer'); cfg.BstSaveTransfer = false; end % 
if ~isfield(cfg,'BstEegTransferFile'); cfg.BstEegTransferFile = 'eeg_transfer.dat'; end % i
if ~isfield(cfg,'BstMegTransferFile'); cfg.BstMegTransferFile = 'meg_transfer.dat'; end % 
if ~isfield(cfg,'BstEegLfFile'); cfg.BstEegLfFile = 'eeg_lf.dat'; end %
if ~isfield(cfg,'BstMegLfFile'); cfg.BstMegLfFile = 'meg_lf.dat'; end % 

