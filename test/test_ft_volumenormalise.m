function test_ft_volumenormalise

% MEM 8gb
% WALLTIME 00:45:00
% DEPENDENCY ft_volumenormalise ft_warp_apply

[ftver, ftpath] = ft_version;

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.mri');
mri = ft_read_mri(filename);

%%
% construct a fake source reconstruction and inperpolate it on the anatomical MRI
% the source reconstruction looks like a checkerboard and is limited to the head

cfg                 = [];
cfg.mri             = mri;
cfg.threshold       = 0.1;
cfg.resolution      = 0.6;
cfg.smooth          = 10;
sourcemodel         = ft_prepare_sourcemodel(cfg);
sourcemodel         = ft_convert_units(sourcemodel, 'mm');

inside = sourcemodel.inside;
dum = zeros(sourcemodel.dim(1:2));
dum(:,1) = 1;
dum(1:2:end,1) = -1;
for k = 2:size(dum,2)
  dum(:,k) = -dum(:,k-1);
end
dum = repmat(dum,[1 1 sourcemodel.dim(3)]);
for k = 2:size(dum,3)
  dum(:,:,k) = -dum(:,:,k-1);
end
dum = dum(:);
dum(~inside) = 0;
sourcemodel.avg.pow = dum;
sourcemodel.inside(:) = true;

cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'nearest';
source = ft_sourceinterpolate(cfg, sourcemodel, mri);

%%

cfg = [];
cfg.funparameter = 'pow';
cfg.location = [0 0 0];
cfg.locationcoordinates = 'head';

cfg.figurename = 'original';
ft_sourceplot(cfg, source);

%%

rmpath(genpath(fullfile(ftpath,'external','spm2'))); rmpath(genpath(fullfile(ftpath,'external','spm8'))); rmpath(genpath(fullfile(ftpath,'external','spm12')));

cfg = [];
cfg.spmversion = 'spm2';
n2 = ft_volumenormalise(cfg, source);

rmpath(genpath(fullfile(ftpath,'external','spm2'))); rmpath(genpath(fullfile(ftpath,'external','spm8'))); rmpath(genpath(fullfile(ftpath,'external','spm12')));

cfg.spmversion = 'spm8';
n8 = ft_volumenormalise(cfg, source);

rmpath(genpath(fullfile(ftpath,'external','spm2'))); rmpath(genpath(fullfile(ftpath,'external','spm8'))); rmpath(genpath(fullfile(ftpath,'external','spm12')));

cfg.spmversion = 'spm12';
cfg.spmmethod = 'old'; % this takes about 11 seconds, it is the default to retain compatibility with older scripts that are based on spm2 or spm8
n12old = ft_volumenormalise(cfg, source);
cfg.spmmethod = 'new'; % this takes about 241 seconds
n12new = ft_volumenormalise(cfg, source);
cfg.spmmethod = 'mars'; % this takes about 210 seconds
n12mars = ft_volumenormalise(cfg, source);

%%

cfg = [];
cfg.funparameter = 'pow';
cfg.location = [0 0 0];
cfg.locationcoordinates = 'head';

cfg.figurename = 'n2';
ft_sourceplot(cfg, n2);
cfg.figurename = 'n8';
ft_sourceplot(cfg, n8);
cfg.figurename = 'n12old';
ft_sourceplot(cfg, n12old);
cfg.figurename = 'n12new';
ft_sourceplot(cfg, n12new);
cfg.figurename = 'n12mars';
ft_sourceplot(cfg, n12mars);

%%
% normalizing the anatomical MRI on the fly should give the same result as using the spmparams in a separate call

rmpath(genpath(fullfile(ftpath,'external','spm2'))); rmpath(genpath(fullfile(ftpath,'external','spm8'))); rmpath(genpath(fullfile(ftpath,'external','spm12')));

cfg = [];
cfg.spmversion = 'spm2';
cfg.spmparams = n2.cfg.spmparams;
n2_params = ft_volumenormalise(cfg, source);
assert(isequal(n2.anatomy, n2_params.anatomy));

rmpath(genpath(fullfile(ftpath,'external','spm2'))); rmpath(genpath(fullfile(ftpath,'external','spm8'))); rmpath(genpath(fullfile(ftpath,'external','spm12')));

cfg = [];
cfg.spmversion = 'spm8';
cfg.spmparams = n8.cfg.spmparams;
n8_params = ft_volumenormalise(cfg, source);
assert(isequal(n8.anatomy, n8_params.anatomy));

rmpath(genpath(fullfile(ftpath,'external','spm2'))); rmpath(genpath(fullfile(ftpath,'external','spm8'))); rmpath(genpath(fullfile(ftpath,'external','spm12')));

cfg = [];
cfg.spmversion = 'spm12';
cfg.spmparams = n12old.cfg.spmparams;
n12old_params = ft_volumenormalise(cfg, source);
assert(isequal(n12old.anatomy, n12old_params.anatomy));

cfg.spmparams = n12new.cfg.spmparams;
n12new_params = ft_volumenormalise(cfg, source);
assert(isequal(n12new.anatomy, n12new_params.anatomy));

cfg.spmparams = n12mars.cfg.spmparams;
n12mars_params = ft_volumenormalise(cfg, source);
assert(isequal(n12mars.anatomy, n12mars_params.anatomy));

% I don't expect the different methods to be exactly the same
assert(~isequal(n12old.anatomy, n12new.anatomy));
assert(~isequal(n12old.anatomy, n12mars.anatomy));


%%
% transition back and forth between individual, approximate and normalised space, using the different sets of parameters

version = {'n2', 'n8', 'n12old', 'n12new', 'n12mars'};

for i=1:numel(version)
  rmpath(genpath(fullfile(ftpath,'external','spm2'))); rmpath(genpath(fullfile(ftpath,'external','spm8'))); rmpath(genpath(fullfile(ftpath,'external','spm12')));
  clear initial params
  
  initial = getfield(eval(version{i}), 'initial');
  params  = getfield(eval(version{i}), 'params');
  
  pos0 = sourcemodel.pos;                                  % individual coordinates
  pos1 = ft_warp_apply(initial, pos0, 'homogenous');       % approximate ACPC coordinates
  pos2 = ft_warp_apply(params, pos1, 'individual2sn');     % MNI coordinates
  pos3 = ft_warp_apply(params, pos2, 'sn2individual');     % approximate ACPC coordinates
  pos4 = ft_warp_apply(inv(initial), pos3, 'homogenous');  % individual coordinates
  
  figure; ft_plot_mesh(pos0(inside, :)); ft_plot_axes([], 'unit', 'mm'); title(sprintf('%s - pos0 - individual coordinates', version{i}))
  figure; ft_plot_mesh(pos1(inside, :)); ft_plot_axes([], 'unit', 'mm'); title(sprintf('%s - pos1 - approximate ACPC coordinates', version{i}))
  figure; ft_plot_mesh(pos2(inside, :)); ft_plot_axes([], 'unit', 'mm'); title(sprintf('%s - pos2 - MNI coordinates', version{i}))
  figure; ft_plot_mesh(pos3(inside, :)); ft_plot_axes([], 'unit', 'mm'); title(sprintf('%s - pos3 - approximate ACPC coordinates', version{i}))
  figure; ft_plot_mesh(pos4(inside, :)); ft_plot_axes([], 'unit', 'mm'); title(sprintf('%s - pos4 - individual coordinates', version{i}))
  
end

%%
% FT_VOLUMENORMALISE is also being called from within FT_PREPARE_SOURCEMODEL

rmpath(genpath(fullfile(ftpath,'external','spm2'))); rmpath(genpath(fullfile(ftpath,'external','spm8'))); rmpath(genpath(fullfile(ftpath,'external','spm12')));

cfg                 = [];
cfg.method          = 'basedonmni';  % this used to be specified with cfg.warpmni = 'yes'
cfg.mri             = mri;
cfg.threshold       = 0.1;
cfg.resolution      = 6;
cfg.nonlinear       = 'yes';
cfg.spmversion      = 'spm12';

cfg.spmmethod = 'old';
sourcemodel_warp_old = ft_prepare_sourcemodel(cfg);

cfg.spmmethod = 'new';
sourcemodel_warp_new = ft_prepare_sourcemodel(cfg);

cfg.spmmethod = 'mars';
sourcemodel_warp_mars = ft_prepare_sourcemodel(cfg);

% they should all be in CTF coordinates, i.e. the frontal pole pointing to the +x axis

ft_determine_coordsys(sourcemodel_warp_old, 'interactive', false)
view(0, 0);

ft_determine_coordsys(sourcemodel_warp_new, 'interactive', false)
view(0, 0);

ft_determine_coordsys(sourcemodel_warp_mars, 'interactive', false)
view(0, 0);
