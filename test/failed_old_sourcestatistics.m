function failed_old_sourcestatistics

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_old_sourcestatistics

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

%this script tests the fixsource function which is part of checkdata

%-------------------------------------
%generate data
rand('twister',20090408);
data = [];
data.fsample = 1000;
data.cfg     = [];

nsmp  = 1000;
nchan = 80;
for k = 1:10
  data.trial{k} = randn(nchan,nsmp);
  data.time{k}  = ([1:nsmp]-1)./data.fsample;
end

%create grad-structure and add to data
[pnt,tri] = icosahedron162;
nrm       = normals(pnt, tri, 'vertex');
pnt       = pnt.*12;
[srt,ind] = sort(pnt(:,3),'descend');
grad      = [];
grad.pnt  = pnt(ind(1:nchan),:);
grad.ori  = nrm(ind(1:nchan),:);
grad.tra  = eye(nchan);
for k = 1:nchan
  grad.label{k} = ['chan',num2str(k,'%03d')];
end
data.grad  = grad;
data.label = grad.label;

%create vol
vol   = [];
vol.o = [0 0 2];
vol.r = 8;

%prepare leadfields and grid
cfg                 = [];
cfg.grid.resolution = 1.5;
cfg.vol             = vol;
cfg.grad            = grad;
grid                = ft_prepare_leadfield(cfg);

%do spectral analysis
cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [2 100];
cfg.pad    = 1;
cfg.tapsmofrq = 3;
freq       = ft_freqanalysis(cfg, data);

%do timelocked analysis
cfg            = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [0 0.999];
cfg.keeptrials = 'yes';
tlck           = ft_timelockanalysis(cfg, data);

%-----------------------------------------------

%-----------------------------------------------
%source reconstruction part

%without fixedori
cfgs           = [];
cfgs.grid      = grid;
cfgs.frequency = 20;
cfgs.method    = 'dics';
cfgs.lambda    = '10%';
cfgs.keepfilter = 'yes';
cfgs.feedback  = 'textbar';
cfgs.vol       = vol;
cfgs.fixedori  = 'yes';
cfgs.keepcsd   = 'yes';
cfgs.realfilter = 'yes';
cfgs.keepleadfield = 'yes';
source         = ft_sourceanalysis(cfgs,freq);

cfgsd            = [];
cfgsd.projectmom = 'yes';
sd               = ft_sourcedescriptives(cfgsd, source);

cfgs.grid.filter = sd.avg.filter;
cfgs.method      = 'pcc';
cfgs.keepmom     = 'yes';

%there are two ways for not crashing the second round of ft_sourceanalysis
%with fixori = 'yes';
if all(islogical(source.inside))
  insidevec = find(source.inside);
else
  insidevec = source.inside;
end

for k = 1:numel(insidevec)
  kk = insidevec(k);
  cfgs.grid.leadfield{kk} = sd.leadfield{kk}*sd.avg.ori{kk};
end
spcc = ft_sourceanalysis(cfgs, freq);

%call sourcedescriptives to get a .trial field
cfgsd            = [];
cfgsd.keeptrials = 'yes';
spcc2            = ft_sourcedescriptives(cfgsd, spcc);

%lcmv
cfgs           = [];
cfgs.grid      = grid;
cfgs.frequency = 20;
cfgs.method    = 'lcmv';
cfgs.lambda    = '10%';
cfgs.keepfilter = 'yes';
cfgs.feedback  = 'textbar';
cfgs.vol       = vol;
cfgs.fixedori  = 'yes';
cfgs.keepcov   = 'yes';
cfgs.keepleadfield = 'yes';
slcmv          = ft_sourceanalysis(cfgs,tlck);

cfgs.grid.filter = slcmv.avg.filter;
cfgs.rawtrial    = 'yes';
slcmv2           = ft_sourceanalysis(cfgs, tlck);

%--------
%do sourcestatistics
cfgs = [];
cfgs.method           = 'montecarlo';
cfgs.parameter        = 'pow';
cfgs.statistic        = 'indepsamplesT';
cfgs.numrandomization = 1;
cfgs.design           = [ones(1,5) ones(1,5)*2];
cfgs.implementation   = 'new';
%statnew1 = ft_sourcestatistics(cfgs, spcc); % this is not going to work,
%because in the new representation, there's no explicit conversion from mom
%into pow
cfgs.implementation   = 'old';
statold1 = ft_sourcestatistics(cfgs, spcc2);
cfgs.parameter        = 'mom';
cfgs.implementation   = 'new';
statnew2 = ft_sourcestatistics(cfgs, slcmv2);
cfgs.implementation   = 'old';
statold2 = ft_sourcestatistics(cfgs, slcmv2);

cfgs = [];
cfgs.method           = 'montecarlo';
cfgs.parameter        = 'pow';
cfgs.statistic        = 'indepsamplesT';
cfgs.numrandomization = 10;
cfgs.correctm         = 'cluster';
cfgs.design           = [ones(1,5) ones(1,5)*2];
cfgs.implementation   = 'new';
cfgs.randomseed       = 20090408;
statnew1 = ft_sourcestatistics(cfgs, spcc);
cfgs.implementation   = 'old';
statold1 = ft_sourcestatistics(cfgs, spcc2);
cfgs.parameter        = 'mom';
cfgs.implementation   = 'new';
statnew2 = ft_sourcestatistics(cfgs, slcmv2);
cfgs.implementation   = 'old';
statold2 = ft_sourcestatistics(cfgs, slcmv2);

%--------
%do sourcestatistics for multiple subjects
sd = ft_sourcedescriptives([], spcc);
sd1 = sd;sd1.avg.pow(sd.inside) = randn(numel(sd.inside),1);
sd2 = sd;sd2.avg.pow(sd.inside) = randn(numel(sd.inside),1);
sd3 = sd;sd3.avg.pow(sd.inside) = randn(numel(sd.inside),1);
sd4 = sd;sd4.avg.pow(sd.inside) = randn(numel(sd.inside),1);
cfgs.parameter      = 'avg.pow';
cfgs.implementation = 'old';
cfgs.statistic      = 'depsamplesT';
cfgs.design         = [1 1 2 2; 1 2 1 2];
cfgs.uvar           = 2;
cfgs.ivar           = 1;
statold3 = ft_sourcestatistics(cfgs, sd1, sd2, sd3, sd4);
cfgs.implementation = 'new';
cfgs.parameter      = 'pow';
statnew3 = ft_sourcestatistics(cfgs, sd1, sd2, sd3, sd4);

%--------
%do sourcestatistics for multiple subjects second level, using descriptive statistic of previous call to sourcestatistics
%FIXME this needs conversion of output of sourcestatistics to old versions


