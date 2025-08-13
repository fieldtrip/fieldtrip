function filename_headmodel = duneuro_write_headmodel(filename, mesh, cfg)
% cfg = bst_prepare_head_model(cfg)
% write the head model file to the current path
% The format can be either ".msh" or ".geo" according
% to the simulation.
% cfg.useTensor = 0;  In this case the conductivity is
%                                 be a scalar, then ".msh" file is used.
% cfg.useTensor = 1;  In this case the conductivity is
%                                 be a tensor, then ".geo" file is used.
%
% File created by Takfarinas MEDANI, November 2019;


% TODO : Optimisation ... use binary input/output
%               Use a standard format for all modalities
%               Discuss with duneuro team to add the I/O files

%% Set the minimal parameters
cfg.useTensor = ft_getopt(cfg, 'useTensor', 0);
cfg.isotrop   = ft_getopt(cfg, 'isotrop',   1);
cfg.layerToKeep = ft_getopt(cfg, 'layerToKeep', []);

if isempty(cfg.layerToKeep)
  cfg.layerToKeep = ones(numel(mesh.tissuelabel), 1);
end

% %% Check if the wole head model will be used in the case of the MEG
% if strcmp(cfg.modality,'meg')
%   % Apply the reduced volume for the MEG, this function will remove the
%   % unselected tissu from the head model, this could be also done for
%   % sEEG later
%   %
%   %%FIXME TO DO
%   cfg = bst_selectVolumeTissu(cfg,cfg.layerToKeep);
% end

if strcmp(cfg.modality,'meeg')
  if sum(cfg.layerToKeep) ~= length(mesh.tissuelabel)
    warning(['THE REDUCED MEG HEAD MODEL WILL BE ALSO USED BY THE EEG'...
      ' >> NOT CORRECT : SEPARATE HEAD MODEL SHOULD BE USED >> NOT IMPLEMENTED YET'])
  end
  %% @@@ WARNING @@@
  % For the combined MEEG .... if the user select specific layers for the MEG
  % the new model will be also used by the EEG .... and this is not correct
  % FOR the MEEG and cfg.layerToKeep contient des zeros .... == > add
  % condition that generate two head model and run duneuro for each
  % modalities. ==> Open an issue ... this could be also a proble for
  % openmeeg
end

if isfield(mesh, 'tet')
  eltype = 'tetrahedron';
elseif isfield(mesh, 'hex')
  eltype = 'hexahedron';
else
  ft_error('unknown element type');
end

%% The file format is dependent on the element type, and on how the conductivity is specified
data = [];
switch eltype
  case 'hexahedron'
    if cfg.isotrop && ~cfg.useTensor
      ext = 'dgf';
      format = 'dgf';
      data = mesh.tissue; % don't subtract 1 on purpose, adjustment for 0-indexing will be done in low-level function out_fem_dgf
    elseif cfg.isotrop && cfg.useTensor
      ft_error('Using the tensor Model with hexahedron mesh is not supported for now');
    elseif ~cfg.isotrop
      ft_error('Using the anisotropy model with hexahedron mesh is not supported for now');
      % TODO : Check if the geo file could be used for the hexa
    end
  case 'tetrahedron'
    if cfg.isotrop && ~cfg.useTensor
      ext = 'msh';
      format = 'gmsh_binary';
      mesh.tetrahedron_regions = mesh.tissue - 1; % subtract 1 on purpose, 0-indexing
    elseif cfg.isotrop && cfg.useTensor
      % Duneuro uses the Cauchy files
      ext = 'geo';
      format = 'geo';
    elseif ~cfg.isotrop && cfg.useTensor
      ext = 'geo';
      format = 'geo';
    elseif ~cfg.isotrop && ~cfg.useTensor
      ft_error('Forbidden combination of options');
    end
end
filename_headmodel = sprintf('%s.%s', filename, ext);
ft_info('writing headmodel geometry to file %s', filename_headmodel);
ft_write_headshape(filename_headmodel, mesh, 'format', format, 'data', data);
