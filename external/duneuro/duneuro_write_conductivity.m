function filename_cond = duneuro_write_conductivity(filename_cond, varargin)

% cfg = bst_prepare_conductivity_model(cfg);
% write the conductivity file either for isotrpic modelor anisotrpic

% Author : Takfarinas MEDANI, October 2019,

% TODO : Optimisation ... use binary input/output
%               Use a standard format for all modalities

conductivity = ft_getopt(varargin, 'conductivity');
conductivity_tensor = ft_getopt(varargin, 'conductivity_tensor');
mesh = ft_getopt(varargin, 'mesh');

if ~isempty(mesh)
  if isfield(mesh, 'tet')
    elem = mesh.tet;
  elseif isfield(mesh, 'hex')
    elem = mesh.hex;
  else
    elem = [];
  end
end

if ~isempty(conductivity)
  [p, f, e] = fileparts(filename_cond);
  if isempty(e)
    filename_cond = sprintf('%s.con', filename_cond);
  end
  write_duneuro_conductivity_file(conductivity, filename_cond);
elseif ~isempty(conductivity_tensor) && ~isempty(elem)
  ft_hastoolbox('brainstorm', 1);

  [p, f, e] = fileparts(filename_cond);
  if isempty(e)
    filename_cond = sprintf('%s.knw', filename_cond);
  end
  
  % Transformation matrix  and tensor mapping on each direction
  CondTensor = zeros(size(elem,1), 6) ;
  for ind = 1:size(elem,1)
    temp0 = reshape(conductivity_tensor(ind,:), 3, []);
    T1 = temp0(:,1:3); % get the 3 eigen vectors
    l  =  diag(temp0(:,4)); % get the eigen value as 3x3
    temp = T1 * l * T1'; % reconstruct the tensors
    CondTensor(ind,:) = [temp(1) temp(5) temp(9) temp(4) temp(8) temp(7)]; % this is the right order       
  end
  % write the tensors
  mesh = [];
  mesh.Elements = elem;
  mesh.Tissue       = mesh.tissue;
  mesh.Tissuelabels = mesh.tissuelabel;
  out_fem_knw(FemMat, CondTensor, filename_cond);
end

function write_duneuro_conductivity_file(conductivity,cond_filename)
% write_duneuro_conductivity_file(conductivity_tensor,cond_filename)
% Create a file with .con extension.
% This file is used by the duneuro application 
% input : conductivity : vector containing the conductivity value of
% each layer.
% Isotropic conductivity only, one value per layer. 

% Author : Takfarinas MEDANI, August 2019,

[filepath,name,ext] = fileparts(cond_filename);
if isempty(ext) || ~strcmp(ext,'.con')
    ext = '.con';
end
cond_filename = fullfile(filepath,[name,ext]);
fid = fopen(cond_filename , 'w');
fprintf(fid, '%d\t', conductivity);
fclose(fid);
