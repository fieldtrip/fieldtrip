function test_bug937

% TEST test_bug937
% TEST ft_prepare_mesh ft_prepare_mesh_new

[pnt, tri] = icosahedron162;

% radiuses and origins are defined in mm
svol(1).o = [0,0,0];
svol(1).r = 10;
svol(1).bnd.pnt = svol(1).r*pnt;
svol(1).bnd.tri = tri;

svol(2).o = [0,0,0];
svol(2).r = 50;
svol(2).bnd.pnt = svol(2).r*pnt;
svol(2).bnd.tri = tri;

svol(3).o = [0,0,0];
svol(3).r = 60;
svol(3).bnd.pnt = svol(3).r*pnt;
svol(3).bnd.tri = tri;

% To generate a volume of 3 concentric spheres (works if number of voxels is odd)
% use this code:
% 
% addpath /home/common/matlab/fieldtrip_private/
% res = 1; % in mm
% for i=3:-1:1
%   tmp2 = zeros(151,151,151);
%   xgrid = -svol(i).r:res:svol(i).r;
%   ygrid = xgrid;
%   zgrid = xgrid;
%   [X, Y, Z]  = ndgrid(xgrid, ygrid, zgrid);
%   pos = [X(:) Y(:) Z(:)];
%   [inside] = bounding_mesh(pos, svol(i).bnd.pnt, svol(i).bnd.tri);
%   l = length(xgrid)
%   c = 76; sel = (l-1)./2; % in voxel
%   tmp = reshape(inside,[l l l]); 
%   tmp2(c-sel:c+sel,c-sel:c+sel,c-sel:c+sel) = tmp;
%   MR{i} = tmp2;
% end
% bkgrnd = zeros(151,151,151);
% bkgrnd = MR{1}+MR{2}+MR{3};

fprintf('Loading a volume with a number N = %d of compartments ... ', numel(svol))
load('test_bug937.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start the different methods 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Both volumetric inputs and meshes can generate more than one mesh at a
% time. This means that the processing of an input with 3 compartments will have to call the
% routine only one time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on mri segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert into a compatible structure 
mri = [];
mri.anatomy = bkgrnd;
mri.sphere1 = bkgrnd==1;
mri.sphere2 = bkgrnd==2;
mri.sphere3 = bkgrnd==3;
mri.transform = eye(4);
mri.dim = size(bkgrnd);
cfg = [];
cfg.tissue      = {'sphere1' 'sphere2' 'sphere3'};
cfg.numvertices = 1000;
bnd = ft_prepare_mesh(cfg,mri);
figure,ft_plot_mesh(bnd(1),'edgecolor','k','facecolor','none')
ft_plot_mesh(bnd(2),'edgecolor','g','facecolor','none')
ft_plot_mesh(bnd(3),'edgecolor','r','facecolor','none')
% figure,ft_plot_mesh(bnd(3),'edgecolor','k','facecolor','w'),camlight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on raw segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see test_bug1646.m
% this input no longer supported

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on headshape, can be a file name, 
%  or a set of points/boundaries (like here)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert in a compatible bnd structure
% set of boundaries
clear hsh
for i=1:3
  hsh(i) = svol(i).bnd;
end
cfg = [];
cfg.headshape   = hsh;
cfg.numvertices = 1000;
bnd = ft_prepare_mesh(cfg);
figure,ft_plot_mesh(bnd(1),'edgecolor','k','facecolor','none')
ft_plot_mesh(bnd(2),'edgecolor','g','facecolor','none')
ft_plot_mesh(bnd(3),'edgecolor','r','facecolor','none')

% set of points
for i=1:3
  cfg = [];
  cfg.headshape   = svol(i).bnd.pnt;
  cfg.numvertices = 1000;
  bnd(i) = ft_prepare_mesh(cfg);
end
figure,ft_plot_mesh(bnd(1),'edgecolor','k','facecolor','none')
ft_plot_mesh(bnd(2),'edgecolor','g','facecolor','none')
ft_plot_mesh(bnd(3),'edgecolor','r','facecolor','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on vol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0 % FIXME (see bug notes)
  vol = svol;
  vol = rmfield(vol,'o');
  vol = rmfield(vol,'r');
  cfg = [];
  cfg.numvertices = 1000;
  bnd = ft_prepare_mesh(cfg,vol);
  figure,ft_plot_mesh(bnd(1),'edgecolor','k','facecolor','none')
  ft_plot_mesh(bnd(2),'edgecolor','g','facecolor','none')
  ft_plot_mesh(bnd(3),'edgecolor','r','facecolor','none')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert in a compatible sph structure
sph = [];
for i=1:3
  sph.r(i)   = svol(i).r;
  sph.o(i,:) = svol(i).o;
end
bnd = ft_prepare_mesh(cfg,sph);
figure,ft_plot_mesh(bnd(1),'edgecolor','k','facecolor','none')
ft_plot_mesh(bnd(2),'edgecolor','g','facecolor','none')
ft_plot_mesh(bnd(3),'edgecolor','r','facecolor','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on interactive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Can't test cfg.interactive='yes' in a test script.
% mri = [];
% mri.anatomy = bkgrnd;
% mri.transform = eye(4);
% mri.dim = size(bkgrnd);
% cfg = [];
% cfg.interactive = 'yes';
% bnd = ft_prepare_mesh(cfg,mri);
