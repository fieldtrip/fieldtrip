function duneuro_write_minifile(cfg, filename)

% DUNEURO_WRITE_MINIFILE creates a text-file with the configuration, 
% to be used in combination with Brainstorm's compiled executable

fid = fopen(filename, 'wt+');

% General settings
fprintf(fid, '__name = %s\n\n', filename);
if strcmp(cfg.solver_type, 'cg')
    fprintf(fid, 'type = %s\n', cfg.type);
end
fprintf(fid, 'element_type = %s\n', cfg.element_type);
fprintf(fid, 'solver_type = %s\n',  cfg.solver_type);
fprintf(fid, 'geometry_adapted = %s\n', bool2str(cfg.geometry_adapted));
fprintf(fid, 'tolerance = %d\n', cfg.tolerance);

% [electrodes]
if strcmp(cfg.modality, 'eeg') || strcmp(cfg.modality, 'meeg')
  fprintf(fid, '[electrodes]\n');
  fprintf(fid, 'filename = %s\n', cfg.filename_elecpos);
  fprintf(fid, 'type = %s\n',     cfg.eeg.type);
  fprintf(fid, 'codims = %s\n', '3');
end
% [meg]
if strcmp(cfg.modality, 'meg') || strcmp(cfg.modality, 'meeg')
  fprintf(fid, '[meg]\n');
  fprintf(fid, 'intorderadd = %d\n',  cfg.meg.intorderadd);
  fprintf(fid, 'type = %s\n',         cfg.meg.type);
  fprintf(fid, 'cache.enable = %s\n', bool2str(cfg.meg.enablecache) );
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
fprintf(fid, 'filename = %s\n', cfg.grid_filename);
% [volume_conductor.tensors]
fprintf(fid, '[volume_conductor.tensors]\n');
fprintf(fid, 'filename = %s\n', cfg.tensors_filename);
% [solver]
fprintf(fid, '[solver]\n');
fprintf(fid, 'solver_type = %s\n',         cfg.solver_type);
fprintf(fid, 'preconditioner_type = %s\n', cfg.solver.preconditioner_type);
if strcmp(cfg.solver_type, 'cg')
    fprintf(fid, 'cg_smoother_type = %s\n', cfg.solver.cg_smoother_type);
end
fprintf(fid, 'intorderadd = %d\n', cfg.solver.intorderadd);
% Discontinuous Galerkin
if strcmp(cfg.solver_type, 'dg')
  fprintf(fid, 'dg_smoother_type = %s\n', cfg.solver.dg_smoother_type);
  fprintf(fid, 'scheme = %s\n',           cfg.solver.dg_scheme);
  fprintf(fid, 'penalty = %d\n',          cfg.solver.dg_penalty);
  fprintf(fid, 'edge_norm_type = %s\n',   cfg.solver.dg_edge_norm_type);
  fprintf(fid, 'weights = %s\n',          bool2str(cfg.solver.dg_weights));
  fprintf(fid, 'reduction = %s\n',        bool2str(cfg.solver.dg_reduction));
end
% [solution]
fprintf(fid, '[solution]\n');
fprintf(fid, 'post_process = %s\n',  bool2str(cfg.post_process)); % true/false
fprintf(fid, 'subtract_mean = %s\n', bool2str(cfg.subtract_mean)); % boolean
% [solution.solver]
fprintf(fid, '[solution.solver]\n');
fprintf(fid, 'reduction = %d\n', cfg.reduction);
% [solution.source_model]
fprintf(fid, '[solution.source_model]\n');
fprintf(fid, 'type = %s\n',              cfg.source_model.type);
fprintf(fid, 'initialization = %s\n',    cfg.source_model.initialization);
fprintf(fid, 'intorderadd = %d\n',       cfg.source_model.intorderadd);
fprintf(fid, 'intorderadd_lb = %d\n',    cfg.source_model.intorderadd_lb);
fprintf(fid, 'numberOfMoments = %d\n',   cfg.source_model.numberOfMoments);
fprintf(fid, 'referenceLength = %d\n',   cfg.source_model.referenceLength);
fprintf(fid, 'weightingExponent = %d\n', cfg.source_model.weightingExponent);
fprintf(fid, 'relaxationFactor = %e\n',  cfg.source_model.relaxationFactor);
fprintf(fid, 'mixedMoments = %s\n',      bool2str(cfg.source_model.mixedMoments));
fprintf(fid, 'restrict = %s\n',          bool2str(cfg.source_model.restrict));
% [brainstorm]
fprintf(fid, '[brainstorm]\n');
fprintf(fid, 'modality = %s\n',                cfg.modality);
fprintf(fid, 'output_folder = %s\n',           cfg.outputpath);
fprintf(fid, 'save_eeg_transfer_file = %s\n',  bool2str(cfg.transfer_save));
fprintf(fid, 'save_meg_transfer_file = %s\n',  bool2str(cfg.transfer_save));
fprintf(fid, 'save_meeg_transfer_file = %s\n', bool2str(cfg.transfer_save));
fprintf(fid, 'eeg_transfer_filename = %s\n',   cfg.transfer_eeg);
fprintf(fid, 'meg_transfer_filename = %s\n',   cfg.transfer_meg);
fprintf(fid, 'eeg_leadfield_filename = %s\n',  cfg.lf_eeg);
fprintf(fid, 'meg_leadfield_filename = %s\n',  cfg.lf_meg);
% Close file
fclose(fid);

try
  system(sprintf('chmod +x %s', filename));
end

function str = bool2str(bool)
if bool
  str = 'true';
else
  str = 'false';
end
