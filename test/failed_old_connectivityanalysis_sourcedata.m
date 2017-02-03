function failed_old_connectivityanalysis_sourcedata

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_old_connectivityanalysis_sourcedata

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

%this script tests the functionality of connectivityanalysis with respect
%to source data

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

%there are two ways for not crashing the second round of sourceanalysis
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
spcc  = ft_sourceanalysis(cfgs, freq);
spcc  = checkdata(spcc, 'sourcerepresentation', 'new');
spccs = source2sparse(spcc);

cfgc = [];
cfgc.method = 'coh';
cfgc.refindx = 1:numel(spccs.inside);
scoh  = ft_connectivityanalysis(cfgc, spccs);
scohf = source2full(scoh); 
