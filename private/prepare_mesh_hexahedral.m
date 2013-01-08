function mesh=prepare_mesh_hexahedral(cfg,mri)

% PREPARE_MESH_HEXAHEDRAL
%
% See also PREPARE_MESH_SEGMENTATION, PREPARE_MESH_MANUAL, PREPARE_MESH_HEADSHAPE
%
% Configuration options for generating a regular 3-D grid
%   cfg.tissue = cell with the names of the compartments that should be
%   meshed
%   cfg.resolution = desired resolution of the mesh (standard = 1)
%
% Copyrights (C) 2012, Johannes Vorwerk
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% get the default options
cfg.tissue      = ft_getopt(cfg, 'tissue');
cfg.resolution  = ft_getopt(cfg, 'resolution');

if isempty(cfg.tissue)
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  fn = fieldnames(mri);
  for i=1:numel(fn),if numel(mri.(fn{i}))==prod(mri.dim), segfield=fn{i};end;end
  cfg.tissue=setdiff(unique(mri.(segfield)(:)),0);
end

if ischar(cfg.tissue)
  % it should either be something like {'brain', 'skull', 'scalp'}, or something like [1 2 3]
  cfg.tissue = {cfg.tissue};
end

if iscell(cfg.tissue)
  % the code below assumes that it is a probabilistic representation
  if any(strcmp(cfg.tissue, 'brain'))
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic', 'hasbrain', 'yes');
  else
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic');
  end
else
  % the code below assumes that it is an indexed representation
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
end

if isempty(cfg.resolution)
  warning('Using standard resolution 1 mm')
  cfg.resolution = 1;
end

% do the mesh extraction
% this has to be adjusted for FEM!!!
if iscell(cfg.tissue)
  % this assumes that it is a probabilistic representation
  % for example {'brain', 'skull', scalp'}
  try
    temp = zeros(size(mri.(cfg.tissue{1})(:)));
    for i=1:numel(cfg.tissue)
      temp = [temp,mri.(cfg.tissue{i})(:)];
    end
    [val,seg] = max(temp,[],2);
    seg = seg - 1;
    seg = reshape(seg,mri.dim);
  catch
    error('Please specify cfg.tissue to correspond to tissue types in the segmented MRI')
  end
  tissue = cfg.tissue;
else
  % this assumes that it is an indexed representation
  % for example [3 2 1]
  seg = zeros(mri.dim);
  tissue = {};
  for i=1:numel(cfg.tissue)
    seg = seg + i*(mri.seg==cfg.tissue(i));
    if isfield(mri, 'seglabel')
      try
        tissue{i} = mri.seglabel{cfg.tissue(i)};
      catch
        error('Please specify cfg.tissue to correspond to (the name or number of) tissue types in the segmented MRI')
      end
    else
      tissue{i} = sprintf('tissue %d', i);
    end
  end
end

% ensure that the segmentation is binary and that there is a single contiguous region
% FIXME is this still needed when it is already binary?
%seg = volumethreshold(seg, 0.5, tissue);

% temporary file names for vgrid call
tname         = tempname;
shfile        = [tname '.sh'];
MRfile        = [tname '_in.v'];
meshfile      = [tname '_out.v'];
materialsfile = [tname '.mtr'];

ft_hastoolbox('simbio', 1);
% write the segmented volume in a Vista format .v file
write_vista_vol(size(seg), seg, MRfile);

% write the materials file (assign tissue values to the elements of the FEM grid)
% see tutorial http://www.rheinahrcampus.de/~medsim/vgrid/manual.html
sb_write_materials(materialsfile, 1:numel(cfg.tissue), tissue, cfg.resolution);

% determin the full path to the vgrid executable
ft_hastoolbox('vgrid', 1);
executable = vgrid; % the helper m-file returns the location of the executable

% write the shell file
efid = fopen(shfile, 'w');
fprintf(efid,'#!/usr/bin/env bash\n');
fprintf(efid,[executable ' -in ' MRfile ' -out ' meshfile ' -min ' num2str(cfg.resolution) ' -max ' num2str(cfg.resolution) ' -elem cube -material ' materialsfile ' -smooth shift -shift 0.30\n']);
fclose(efid);
dos(sprintf('chmod +x %s', shfile));
disp('vgrid is writing the wireframe mesh file, this may take some time ...')
stopwatch = tic;

% use vgrid to construct the wireframe
system(shfile);
fprintf('elapsed time: %d seconds\n', toc(stopwatch));


% read the mesh points
[mesh.pnt,mesh.hex,labels] = read_vista_mesh(meshfile);

% converting position of meshpoints to the head coordinate system

mesh.pnt = warp_apply(mri.transform,mesh.pnt,'homogeneous');


mesh.tissue = zeros(size(labels));
numlabels = size(unique(labels),1);
mesh.tissuelabel = {};
ulabel = sort(unique(labels));
for i = 1:numlabels
  mesh.tissue(labels == ulabel(i)) = i;
  mesh.tissuelabel{i} = num2str(ulabel(i));
  mesh.tissuename{i} = tissue{i};
end

end % function
