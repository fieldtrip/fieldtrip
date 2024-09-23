function duneuro_write_minifile(cfg, IniFile)

%% ===== WRITE MINI FILE =====
% Open the mini file
fid = fopen(IniFile, 'wt+');

% General settings
fprintf(fid, '__name = %s\n\n', IniFile);
if strcmp(cfg.dnFemSolverType, 'cg')
    fprintf(fid, 'type = %s\n', cfg.dnFemMethodType);
end
fprintf(fid, 'element_type = %s\n', cfg.dnMeshElementType);
fprintf(fid, 'solver_type = %s\n', cfg.dnFemSolverType);
fprintf(fid, 'geometry_adapted = %s\n', bool2str(cfg.dnGeometryAdapted));
fprintf(fid, 'tolerance = %d\n', cfg.dnTolerance);

% [electrodes]

%%FIXME
isEcog = false; isSeeg = false;
if isEcog || isSeeg
  % Instead of selecting the electrode on the outer surface,
  % uses the nearest FEM node as the electrode location
  cfg.ElecType = 'closest_subentity_center';
end
if strcmp(cfg.modality, 'eeg') || strcmp(cfg.modality, 'meeg')
  fprintf(fid, '[electrodes]\n');
  fprintf(fid, 'filename = %s\n', cfg.filename_elecpos);
  fprintf(fid, 'type = %s\n',     cfg.dnElectrodType);
  fprintf(fid, 'codims = %s\n', '3');
end
% [meg]
if strcmp(cfg.modality, 'meg') || strcmp(cfg.modality, 'meeg')
  fprintf(fid, '[meg]\n');
  fprintf(fid, 'intorderadd = %d\n',  cfg.dnMegIntorderadd);
  fprintf(fid, 'type = %s\n',         cfg.dnMegType);
  fprintf(fid, 'cache.enable = %s\n', bool2str(cfg.EnableCacheMemory) );
  % [coils]
  fprintf(fid, '[coils]\n');
  fprintf(fid, 'filename = %s\n', cfg.filename_coilpos);
  % [projections]
  fprintf(fid, '[projections]\n');
  fprintf(fid, 'filename = %s\n', cfg.filename_coilori);
end
% [dipoles]
fprintf(fid, '[dipoles]\n');
fprintf(fid, 'filename = %s\n', cfg.filename_dipoles);
% [volume_conductor.grid]
fprintf(fid, '[volume_conductor.grid]\n');
fprintf(fid, 'filename = %s\n', cfg.filename_headmodel);
% [volume_conductor.tensors]
fprintf(fid, '[volume_conductor.tensors]\n');
fprintf(fid, 'filename = %s\n', cfg.filename_cond);
% [solver]
fprintf(fid, '[solver]\n');
fprintf(fid, 'solver_type = %s\n',         cfg.dnSolverSolverType);
fprintf(fid, 'preconditioner_type = %s\n', cfg.dnSolverPreconditionerType);
if strcmp(cfg.dnSolverSolverType, 'cg')
    fprintf(fid, 'cg_smoother_type = %s\n', cfg.dnSolverCgSmootherType);
end
fprintf(fid, 'intorderadd = %d\n', cfg.dnSolverIntorderadd);
% Discontinuous Galerkin
if strcmp(cfg.dnSolverSolverType, 'dg')
  fprintf(fid, 'dg_smoother_type = %s\n', cfg.DgSmootherType);
  fprintf(fid, 'scheme = %s\n',           cfg.DgScheme);
  fprintf(fid, 'penalty = %d\n',          cfg.DgPenalty);
  fprintf(fid, 'edge_norm_type = %s\n',   cfg.DgEdgeNormType);
  fprintf(fid, 'weights = %s\n',          bool2str(cfg.DgWeights));
  fprintf(fid, 'reduction = %s\n',        bool2str(cfg.DgReduction));
end
% [solution]
fprintf(fid, '[solution]\n');
fprintf(fid, 'post_process = %s\n', bool2str(cfg.dnSolutionPostProcess)); % true/false
fprintf(fid, 'subtract_mean = %s\n', bool2str(cfg.dnSolutionSubstractMean)); % boolean
% [solution.solver]
fprintf(fid, '[solution.solver]\n');
fprintf(fid, 'reduction = %d\n', cfg.dnSolutionSolverReduction);
% [solution.source_model]
fprintf(fid, '[solution.source_model]\n');
fprintf(fid, 'type = %s\n',              cfg.femSourceModel);
fprintf(fid, 'intorderadd = %d\n',       cfg.femSourceModelIntorderadd);
fprintf(fid, 'intorderadd_lb = %d\n',    cfg.femSourceModelIntorderadd_lb);
fprintf(fid, 'numberOfMoments = %d\n',   cfg.femSourceModelNumberOfMoments);
fprintf(fid, 'referenceLength = %d\n',   cfg.femSourceModelReferenceLength);
fprintf(fid, 'weightingExponent = %d\n', cfg.femSourceModelWeightingExponent);
fprintf(fid, 'relaxationFactor = %e\n',  10^(-cfg.femSourceModelRelaxationFactor));
fprintf(fid, 'mixedMoments = %s\n',      bool2str(cfg.femSourceModelMixedMoments));
fprintf(fid, 'restrict = %s\n',          bool2str(cfg.femSourceModelRestrict));
fprintf(fid, 'initialization = %s\n',    cfg.femSourceModelInitialization);
% [brainstorm]
fprintf(fid, '[brainstorm]\n');
fprintf(fid, 'modality = %s\n',                cfg.modality);
fprintf(fid, 'output_folder = %s\n',           cfg.BstOutputFolder);
fprintf(fid, 'save_eeg_transfer_file = %s\n',  bool2str(cfg.BstSaveTransfer));
fprintf(fid, 'save_meg_transfer_file = %s\n',  bool2str(cfg.BstSaveTransfer));
fprintf(fid, 'save_meeg_transfer_file = %s\n', bool2str(cfg.BstSaveTransfer));
fprintf(fid, 'eeg_transfer_filename = %s\n',   cfg.BstEegTransferFile);
fprintf(fid, 'meg_transfer_filename = %s\n',   cfg.BstMegTransferFile);
fprintf(fid, 'eeg_leadfield_filename = %s\n',  cfg.BstEegLfFile);
fprintf(fid, 'meg_leadfield_filename = %s\n',  cfg.BstMegLfFile);
% Close file
fclose(fid);

try
  system(sprintf('chmod +x %s', IniFile));
end

function str = bool2str(bool)
if bool
  str = 'true';
else
  str = 'false';
end
