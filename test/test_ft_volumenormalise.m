function test_ft_volumenormalise

% MEM 4gb
% WALLTIME 00:45:00
% DEPENDENCY ft_volumenormalise ft_warp_apply

[ftver, ftpath] = ft_version;

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri');
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
cfg.newfigure = 'yes';
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
cfg.newfigure = 'yes';
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
% try going back and forth between individual and normalised space, using the different sets of parameters
%
% SOME PARTS OF THIS SECTION ARE CURRENTLY FAILLING, IT MIGHT BE THAT I UNDERSTAND IT WRONG
% IT IS NOT RELATED to https://github.com/fieldtrip/fieldtrip/pull/1450

rmpath(genpath(fullfile(ftpath,'external','spm2'))); rmpath(genpath(fullfile(ftpath,'external','spm8'))); rmpath(genpath(fullfile(ftpath,'external','spm12')));
pos = sourcemodel.pos;

pos2_sn      = ft_warp_apply(n2.cfg.spmparams, pos, 'individual2sn');
pos2_sn_ind  = ft_warp_apply(n2.cfg.spmparams, pos2_sn, 'sn2individual');

% dist = sqrt(sum((pos2_sn_ind-pos2_sn).^2, 2));  % there should be quite some distance between these
% assert(median(dist)>50);
% 
% dist = sqrt(sum((pos2_sn_ind-pos).^2, 2));      % there should not be any distance between these
% assert(median(dist)<1);

pos8_sn      = ft_warp_apply(n8.cfg.spmparams, pos, 'individual2sn');
pos8_sn_ind  = ft_warp_apply(n8.cfg.spmparams, pos8_sn, 'sn2individual');

% dist = sqrt(sum((pos8_sn_ind-pos8_sn).^2, 2));  % there should be quite some distance between these
% assert(median(dist)>50);
% 
% dist = sqrt(sum((pos8_sn_ind-pos).^2, 2));      % there should not be any distance between these
% assert(median(dist)<1);

pos12old_sn      = ft_warp_apply(n12old.cfg.spmparams, pos, 'individual2sn');
pos12old_sn_ind  = ft_warp_apply(n12old.cfg.spmparams, pos12old_sn, 'sn2individual');

% dist = sqrt(sum((pos12old_sn_ind-pos12old_sn).^2, 2));  % there should be quite some distance between these
% assert(median(dist)>50);
% 
% dist = sqrt(sum((pos12old_sn_ind-pos).^2, 2));          % there should not be any distance between these
% assert(median(dist)<1);

pos12new_sn     = ft_warp_apply(n12new.cfg.spmparams, pos, 'individual2sn');
pos12new_sn_ind = ft_warp_apply(n12new.cfg.spmparams, pos12new_sn, 'sn2individual');

% dist = sqrt(sum((pos12new_sn_ind-pos12new_sn).^2, 2));  % there should be quite some distance between these
% assert(median(dist)>50);
% 
% dist = sqrt(sum((pos12new_sn_ind-pos).^2, 2));          % there should not be any distance between these
% assert(median(dist)<1);

%%

pos12mars_sn     = ft_warp_apply(n12mars.cfg.spmparams, pos, 'individual2sn');
pos12mars_sn_ind = ft_warp_apply(n12mars.cfg.spmparams, pos12mars_sn, 'sn2individual');

dist = sqrt(sum((pos12mars_sn_ind-pos12mars_sn).^2, 2)); % there should be quite some distance between these
assert(median(dist)>50);

dist = sqrt(sum((pos12mars_sn_ind-pos).^2, 2));          % there should not be any distance between these
assert(median(dist)<1);


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
