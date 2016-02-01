function test_bug1227

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1227 ft_sourcewrite

[pnt, tri] = icosahedron162;

source = [];
source.pos = pnt;
source.tri = tri;
source.pow = randn(162,100);
source.dimord = 'pos_time';
source.time = (0:99)/100;

cfg = [];
cfg.filename  = [tempname,'.gii'];
cfg.parameter = 'pow'; 
ft_sourcewrite(cfg, source);
if exist(cfg.filename, 'file')
  fprintf('gifti file created\n');
  mrig = ft_read_headshape(cfg.filename);
  delete(cfg.filename);
  fprintf('gifti file deleted\n');
end

[x, y, z]    = ndgrid(1:4,1:5,1:6);
source2      = [];
source2.pos  = [x(:) y(:) z(:)];
source2.dim  = [4 5 6];
source2.time = (0:99)/100;
source2.pow  = randn(120,100);

cfg.filename = [tempname,'.nii'];
ft_sourcewrite(cfg, source2);
if exist(cfg.filename, 'file')
  fprintf('nifti file created\n');
  mrin = ft_read_mri(cfg.filename);
  delete(cfg.filename);
  fprintf('nifti file deleted\n');
end
