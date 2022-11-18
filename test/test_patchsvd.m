function test_patchsvd(datadirs)

% WALLTIME 00:30:00
% MEM 6gb
% DEPENDENCY ft_prepare_leadfield

%if nargin==0
    datadirs{1} = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf');
    datadirs{2} = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/networkanalysis');
%end

%% read the continuous data and segment into 2 seconds epochs, with 50% overlap
cfg            = [];
cfg.dataset    = fullfile(datadirs{1},'SubjectRest.ds'); 
cfg.continuous = 'yes';
cfg.channel    = {'MEG'};
data1 = ft_preprocessing(cfg);

cfg         = [];
cfg.length  = 2;
cfg.overlap = 0;
data1       = ft_redefinetrial(cfg, data1);

data1 = removefields(data1, 'elec');

datadir = datadirs{2};
cd(datadir);


%% downsample the data to speed up component analysis
data1.time(1:end) = data1.time(1);

cfg            = [];
cfg.resamplefs = 100;
cfg.detrend    = 'yes';
data1         = ft_resampledata(cfg, data1);

%% load the required geometrical information
load(fullfile(datadir, 'hdm.mat'));

data1.grad = ft_convert_units(data1.grad, 'm');
hdm = ft_convert_units(hdm, 'm');

%% compute sourcemodels
cfg            = [];
cfg.method     = 'basedonresolution';
cfg.resolution = 0.004;
cfg.headmodel  = hdm;
sourcemodel4   = ft_prepare_sourcemodel(cfg);
%cfg.resolution = 0.008;
%sourcemodel8   = ft_prepare_sourcemodel(cfg);

sel = false(sourcemodel4.dim);
sel(1:2:end,1:2:end,1:2:end) = true;
sourcemodel8 = sourcemodel4;
sourcemodel8.pos = sourcemodel8.pos(sel(:),:);
sourcemodel8.inside = sourcemodel8.inside(sel(:));
sourcemodel8.dim = ceil(sourcemodel4.dim/2);

atlas = sourcemodel8;

p = zeros(atlas.dim);
p = repmat((1:atlas.dim(1))',[1 atlas.dim(2:end)]);
p(~atlas.inside) = 0;

atlas.parcellation = p(:);
for k = 1:atlas.dim(1)
  atlas.parcellationlabel{k,1} = sprintf('parcel%02d',k);
end

%% compute the leadfield
cfg             = [];
cfg.sourcemodel = sourcemodel8;
cfg.headmodel   = hdm;
cfg.channel     = {'MEG'};
cfg.singleshell.batchsize = 1500;
lf            = ft_prepare_leadfield(cfg, data1);

cfg.patchsvd = 0.03;
cfg.patchsvdnum = 5;
lfpatch = ft_prepare_leadfield(cfg, data1);

cfg = removefields(cfg, {'patchsvd' 'patchsvdnum'});
cfg.patchsvd = 'yes';
cfg.atlas    = atlas;
cfg.parcellation = 'parcellation';
lfpatch2 = ft_prepare_leadfield(cfg, data1);
