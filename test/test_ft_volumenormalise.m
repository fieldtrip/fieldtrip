function test_ft_volumenormalise

% MEM 4gb
% WALLTIME 00:45:00

% DEPENDENCY ft_volumenormalise ft_warp_apply

filename = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri');
mri = ft_read_mri(filename);

% try going back and forth (between individual and normalised space), using
% the different sets of parameters
cfg                 = [];
cfg.mri             = mri;
cfg.threshold       = 0.1;
cfg.sourcemodel.resolution = 0.6;
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

ft_file = which('ft_defaults');
[p,f,e] = fileparts(ft_file);

cfg = [];
cfg.spmversion = 'spm2';
n2 = ft_volumenormalise(cfg, source);

rmpath(fullfile(p,'external','spm2'));
cfg.spmversion = 'spm8';
n8 = ft_volumenormalise(cfg, source);

rmpath(fullfile(p,'external','spm8'));
cfg.spmversion = 'spm12';
n12 = ft_volumenormalise(cfg, source);

cfg.spmmethod = 'new';
n12new = ft_volumenormalise(cfg, source);


rmpath(genpath(fullfile(p,'external','spm12')));
P        = n12.cfg.spmparams;
pos      = sourcemodel.pos;
pos12_sn = ft_warp_apply(P, pos, 'individual2sn');
pos12_sn_ind = ft_warp_apply(P, pos12_sn, 'sn2individual');



P        = n12new.cfg.spmparams;
pos      = sourcemodel.pos;
pos12new_sn = ft_warp_apply(P, pos, 'individual2sn');
pos12new_sn_ind = ft_warp_apply(P, pos12new_sn, 'sn2individual');


cfg                 = [];
cfg.spmversion      = 'spm12';
cfg.spmmethod       = 'new';
cfg.mri             = mri;
cfg.threshold       = 0.1;
cfg.sourcemodel.resolution = 6;
cfg.sourcemodel.warpmni    = 'yes';
cfg.sourcemodel.nonlinear  = 'yes';
sourcemodel_warp    = ft_prepare_sourcemodel(cfg);

rmpath(genpath(fullfile(p,'external','spm12')));
cfg.spmmethod = 'old';
sourcemodel_warp2 = ft_prepare_sourcemodel(cfg);
