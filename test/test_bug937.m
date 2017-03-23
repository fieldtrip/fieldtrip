function test_bug937

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_mesh ft_voltype ft_plot_mesh


csvol.o = [0,0,0];
csvol.r = [10 50 60];
assert(ft_voltype(csvol,'concentricspheres'))

ssvol.o = [0,0,0];
ssvol.r = [60];
assert(ft_voltype(ssvol,'singlesphere'))

[pnt, tri] = icosahedron162;
svol.bnd(1).pnt = 10*pnt;
svol.bnd(1).tri = tri;
svol.bnd(2).pnt = 50*pnt;
svol.bnd(2).tri = tri;
svol.bnd(3).pnt = 60*pnt;
svol.bnd(3).tri = tri;
assert(ft_voltype(svol,'unknown'))

tcfg=[];
tcfg.headshape=svol.bnd;
svolcs=ft_prepare_concentricspheres(tcfg);
assert(ft_voltype(svolcs,'concentricspheres'))


% To generate a volume of 3 concentric spheres (works if number of voxels is odd)
% use this code:
%
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

% fprintf('Loading a volume with a number N = %d of compartments ... ', numel(svol))
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug937.mat'));

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
  hsh(i) = svol.bnd(i);
end
cfg = [];
cfg.headshape   = hsh;
% cfg.headshape   = svol;
cfg.numvertices = 1000;
bnd = ft_prepare_mesh(cfg);
figure,ft_plot_mesh(bnd(1),'edgecolor','k','facecolor','none')
ft_plot_mesh(bnd(2),'edgecolor','g','facecolor','none')
ft_plot_mesh(bnd(3),'edgecolor','r','facecolor','none')

% set of points
clear bnd
for i=1:3
  cfg = [];
  cfg.headshape   = svol.bnd(i).pnt;
  cfg.numvertices = 1000;
  bnd(i) = ft_prepare_mesh(cfg);
end
figure,ft_plot_mesh(bnd(1),'edgecolor','k','facecolor','none')
ft_plot_mesh(bnd(2),'edgecolor','g','facecolor','none')
ft_plot_mesh(bnd(3),'edgecolor','r','facecolor','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on vol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
  cfg = [];
  cfg.numvertices = 1000;
  bnd = ft_prepare_mesh(cfg,svol);
  success=true;
catch
  success=false;
end
if success,error('this svol should not work here');end

cfg = [];
cfg.numvertices = 1000;
% convert in a compatible sphere structure
bnd = ft_prepare_mesh(cfg,csvol);
figure,ft_plot_mesh(bnd(1),'edgecolor','k','facecolor','none')
ft_plot_mesh(bnd(2),'edgecolor','g','facecolor','none')
ft_plot_mesh(bnd(3),'edgecolor','r','facecolor','none')
bnd = ft_prepare_mesh(cfg,svolcs);
figure,ft_plot_mesh(bnd(1),'edgecolor','k','facecolor','none')
ft_plot_mesh(bnd(2),'edgecolor','g','facecolor','none')
ft_plot_mesh(bnd(3),'edgecolor','r','facecolor','none')
% convert in a compatible sphere structure
bnd = ft_prepare_mesh(cfg,ssvol);
figure,ft_plot_mesh(bnd,'edgecolor','k','facecolor','none')


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
