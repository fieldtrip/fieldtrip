function cfg = duneuro_defaults(cfg)

% This function is based on the bst-duneuro functon bst_load_default_duneuroConfiguration
 
if nargin<1
  cfg = [];
end

%% General settings
cfg.type                           = ft_getopt(cfg, 'type',         'fitted'); % 'fitted' or 'unfitted'
cfg.solver_type                    = ft_getopt(cfg, 'solver_type',  'cg');     % what else?
cfg.element_type                   = ft_getopt(cfg, 'element_type', []);       % this should be determined from the mesh
cfg.application                    = ft_getopt(cfg, 'application',  []);

cfg.dnGeometryAdapted              = ft_getopt(cfg, 'dnGeometryAdapted',              false);  % true or  false why and  how
cfg.dnTolerance                    = ft_getopt(cfg, 'dnTolerance',                    1e-8);

%% 1 Sensors
% subpart electrode : [electrodes]
cfg.eeg               = ft_getopt(cfg,     'eeg');
cfg.eeg.type          = ft_getopt(cfg.eeg, 'type',  'closest_subentity_center');
cfg.eeg.subentities   = ft_getopt(cfg.eeg, 'subentities', [1 2 3]);

%% subpart [meg]
cfg.meg               = ft_getopt(cfg,     'meg');
cfg.meg.intorderadd   = ft_getopt(cfg.meg, 'intorderadd',   2); % FIXME (was 0 in bst, in FT seems to have 2 as default): why is this numeric, is it intended boolean? why also false/true as strings?
cfg.meg.type          = ft_getopt(cfg.meg, 'type',          'physical');
cfg.meg.enablecache   = ft_getopt(cfg.meg, 'enablecache',   false); % no idea

%% 4 - Subpart  [solver] ==> refers to the linear system solver ?
cfg.solver                     = ft_getopt(cfg,        'solver');
cfg.solver.preconditioner_type = ft_getopt(cfg.solver, 'preconditioner_type', 'amg');  % what are the others 
cfg.solver.cg_smoother_type    = ft_getopt(cfg.solver, 'cg_smoother_type',    'ssor'); % what are the others 
cfg.solver.intorderadd         = ft_getopt(cfg.solver, 'intorderadd', 2);

% case of the dg discontinious galerkin
cfg.solver.dg_smoother_type  = ft_getopt(cfg.solver, 'dg_smoother_type', 'ssor');
cfg.solver.dg_scheme         = ft_getopt(cfg.solver, 'dg_scheme', 'sipg');
cfg.solver.dg_penalty        = ft_getopt(cfg.solver, 'dg_penalty', 20);
cfg.solver.dg_edge_norm_type = ft_getopt(cfg.solver, 'dg_edge_norm_type', 'houston');
cfg.solver.dg_weights        = ft_getopt(cfg.solver, 'dg_weights',   true);
cfg.solver.dg_reduction      = ft_getopt(cfg.solver, 'dg_reduction', true);

%% 5 - Subpart  [solution]
cfg.post_process  = ft_getopt(cfg, 'post_process',  true);
cfg.subtract_mean = ft_getopt(cfg, 'subtract_mean', true); 

% subpart  [solution.solver]
cfg.reduction    = ft_getopt(cfg, 'reduction', 1e-15);

%% 6 - subpart  [solution.source_model]
cfg.source_model                  = ft_getopt(cfg,              'source_model');
cfg.source_model.type             = ft_getopt(cfg.source_model, 'type', 'venant'); % partial_integration, venant, subtraction | expand smtype
cfg.source_model.initialization   = ft_getopt(cfg.source_model, 'initialization', 'closest_vertex');
cfg.source_model.intorderadd      = ft_getopt(cfg.source_model, 'intorderadd', 2);
cfg.source_model.intorderadd_lb   = ft_getopt(cfg.source_model, 'intorderadd_lb', 3);
cfg.source_model.numberOfMoments  = ft_getopt(cfg.source_model, 'numberOfMoments', 3); 
cfg.source_model.referenceLength  = ft_getopt(cfg.source_model, 'referenceLength', 20); 
cfg.source_model.relaxationFactor = ft_getopt(cfg.source_model, 'relaxationFactor', 1e-6); 
cfg.source_model.restrict         = ft_getopt(cfg.source_model, 'restrict', true); 
cfg.source_model.weightingExponent = ft_getopt(cfg.source_model, 'weightingExponent', 1); 
cfg.source_model.mixedMoments     = ft_getopt(cfg.source_model, 'mixedMoments', true); 

% figure out whether the user provided an application explicitly
if isempty(cfg.application) && ~isfield(cfg, 'bstflag')
  % default to duneuro_meeg, which relies on a platform-specific compiled mex-file duneuro_matlab.mex**
  cfg.bstflag = false;
elseif exist(cfg.application, 'file') && contains(cfg.application, 'bst_duneuro')
  % assume it to be a valid compiled application downloaded from brainstorm
  cfg.bstflag = true;
end

if cfg.bstflag
  cfg.outputpath = ft_getopt(cfg, 'outputpath',  tempname);
  if ~endsWith(cfg.outputpath, filesep)
    cfg.outputpath = sprintf('%s%s', cfg.outputpath,filesep); % somehow this is needed 
  end
  if ~exist(cfg.outputpath, 'dir')
    mkdir(cfg.outputpath);
  end
  cfg.minifile_filename = ft_getopt(cfg, 'minifile_filename',  fullfile(cfg.outputpath, 'duneuro_minifile.mini'));

  cfg.modality = ft_getopt(cfg, 'modality', []); % this cannot be guessed and needs to be figured out elsewhere

  %% 7 - subpart  [brainstorm], specify the filenames and some other things
  cfg.transfer_save = ft_getopt(cfg, 'transfer_save', false);
  cfg.transfer_eeg  = ft_getopt(cfg, 'transfer_eeg',  'eeg_transfer.dat');
  cfg.transfer_meg  = ft_getopt(cfg, 'transfer_meg',  'meg_transfer.dat');
  cfg.lf_eeg        = ft_getopt(cfg, 'lf_eeg', 'eeg_lf.dat');
  cfg.lf_meg        = ft_getopt(cfg, 'lf_meg', 'meg_lf.dat');
end
