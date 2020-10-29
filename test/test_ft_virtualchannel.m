function test_ft_virtualchannel

% WALLTIME 00:45:00
% MEM 6gb
% DEPENDENCY ft_virtualchannel

%% do some virtual channel stuff
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/source_coh_lft.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/source_diff.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/data_cmb.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/sourcemodel.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/hdm.mat'));

[maxval, maxcohindx] = max(source_coh_lft.avg.coh);
[maxval, maxpowindx] = max(source_diff.avg.pow);

cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = 'MEG';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
tlock                 = ft_timelockanalysis(cfg, data_cmb);

cfg        = [];
cfg.method = 'mtmfft';
cfg.foi    = 10;
cfg.taper  = 'hanning';
cfg.output = 'fourier';
freq       = ft_freqanalysis(cfg, data_cmb);

% this is old-style stuff. as of end 2020 there's a ft_virtualchannel
% function that does the virtualchannel creation
cfg             = [];
cfg.headmodel   = hdm;
cfg.sourcemodel = sourcemodel;
cfg.channel     = 'MEG';
cfg.singleshell.batchsize = 2000;
leadfield       = ft_prepare_leadfield(cfg, tlock);

cfg              = [];
cfg.method       = 'lcmv';
cfg.sourcemodel  = leadfield;
cfg.lcmv.keepfilter = 'yes';
source           = ft_sourceanalysis(cfg, tlock);

%% create virtual channel data using different methods, and do a few sanity checks
cfg                = [];
cfg.pos            = source.pos([maxcohindx maxpowindx],:);
cfg.method         = 'svd';
cfg.numcomponent   = 1;
data_vc1           = ft_virtualchannel(cfg, data_cmb, source);
data_vc1f          = ft_virtualchannel(cfg, freq, source);

cfg                = [];
cfg.pos            = source.pos([maxcohindx maxpowindx],:);
cfg.method         = 'pca';
cfg.numcomponent   = 1;
data_vc2           = ft_virtualchannel(cfg, data_cmb, source);
data_vc2b          = ft_virtualchannel(cfg, tlock, source);

% data_vc1 and data_vc2 should be very similar, up to a polarity difference
% this can be checked from the balancing, as well as from the time domain
% data
dat1 = cat(2, data_vc1.trial{:});
dat2 = cat(2, data_vc2.trial{:});
cmat = corr([dat1' dat2']);
cmat = cmat([1 3],[1 3]);
assert(all(abs(cmat(:))>0.999));

tra1 = data_vc1.grad.balance.virtualchannel.tra(1:2,:);
tra2 = data_vc2.grad.balance.virtualchannel.tra(1:2,:);
cmat = corr([tra1' tra2']);
cmat = cmat([1 3],[1 3]);
assert(all(abs(cmat(:))>0.99));

% data_vc1f should be of type freq
assert(ft_datatype(data_vc1f, 'freq'));

% data_vc2b should be of type timelock
assert(ft_datatype(data_vc2b, 'timelock'));

cfg                = [];
cfg.pos            = source.pos([maxcohindx maxpowindx],:);
cfg.method         = 'svd';
cfg.numcomponent   = 3;
data_vc3           = ft_virtualchannel(cfg, data_cmb, source);

% data_vc3 should have each vc's third row to be ~0 due to the rank 2
% spatial filters
dat3 = cat(2, data_vc3.trial{:});
vdat3 = var(dat3, [], 2);
vdat3(1:3,:) = vdat3(1:3)./vdat3(1);
vdat3(4:6,:) = vdat3(4:6)./vdat3(4);
assert(all(vdat3([3 6])<1e-12));

cfg                = [];
cfg.pos            = source.pos([maxcohindx maxpowindx],:);
cfg.method         = 'none';
data_vc4           = ft_virtualchannel(cfg, data_cmb, source);

% data_vc4 should have a balancing array that is the same as the spatial
% filters (up to the order of the channels)
f1 = source.avg.filter{maxcohindx};
f4 = full(data_vc4.grad.balance.virtualchannel.tra(1:3,1:151));
assert(isequal(f1,f4));

% ensure that a random order of the channels does not affect the output
ix = randperm(numel(data_cmb.label));
data_cmb2 = data_cmb;
data_cmb2.label = data_cmb2.label(ix);
for k = 1:numel(data_cmb2.trial)
  data_cmb2.trial{k} = data_cmb2.trial{k}(ix,:);
end
cfg                = [];
cfg.pos            = source.pos([maxcohindx maxpowindx],:);
cfg.method         = 'svd';
cfg.numcomponent   = 1;
data_vc1b          = ft_virtualchannel(cfg, data_cmb2, source);

dat1 = cat(2, data_vc1.trial{:});
dat1b = cat(2, data_vc1b.trial{:});
cmat = corr([dat1' dat1b']);
cmat = cmat([1 3],[1 3]);
assert(all(abs(cmat(:))>0.999));

% try a simple parcellation, consisting of 2 parcels, 
tmp = load('standard_sourcemodel3d8mm.mat');
s   = tmp.sourcemodel;
s.parcellation = zeros(size(s.pos,1),1);
s.parcellation(maxpowindx) = 1;
s.parcellation(maxcohindx) = -1;

cur_dir = pwd;
[a,b] = ft_version;
cd(fullfile(b, 'private'));
s.parcellation = reshape(volumesmooth(reshape(s.parcellation,s.dim), 2), [], 1);
cd(cur_dir);

sel1 = s.parcellation>0.02;
sel2 = s.parcellation<-0.02;
s.parcellation(:) = 0;
s.parcellation(sel1) = 1;
s.parcellation(sel2) = 2;
s.parcellationlabel  = {'visual'; 'motor'};
s.pos = source.pos;

cfg = [];
cfg.method = 'svd';
cfg.numcomponent = 1;
cfg.parcellation = 'parcellation';
data_vc1p = ft_virtualchannel(cfg, data_cmb, source, s);

% compute power constrast vs. bsl.
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [0 100];
cfg.pad    = 1;
cfg.tapsmofrq = 5;
freq = ft_freqanalysis(cfg, data_vc1p);

cfg = [];
cfg.trials = find(freq.trialinfo == 0);
bsl = ft_freqdescriptives(cfg, freq);
cfg.trials = find(freq.trialinfo == 1);
act = ft_freqdescriptives(cfg, freq);

figure; plot(bsl.freq, act.powspctrm./bsl.powspctrm - 1);

ft_hastoolbox('dss', 1);
cfg               = [];
cfg.cellmode      = 'yes';
cfg.method        = 'dss';
cfg.demean        = 'no';
cfg.doscale       = 'yes';
cfg.dss.algorithm = 'pca';
cfg.dss.denf.function = 'denoise_filter2';
cfg.dss.denf.params.filter_filtfilt.A = [];
cfg.dss.denf.params.filter_filtfilt.B = [];

[dum, B, A] = ft_preproc_bandpassfilter(data_cmb.trial{end},data_cmb.fsample,[40 60],[],'firws');
cfg.dss.denf.params.filter_filtfilt.A = A;
cfg.dss.denf.params.filter_filtfilt.B = B;
cfg.dss.denf.params.filter_filtfilt.function = 'fir_filterdcpadded';

%cfg.trials = find(data_cmb.trialinfo==1);
cfg.method = 'dss';
cfg.parcellation = 'parcellation';
data_vc2p = ft_virtualchannel(cfg, data_cmb, source, s);

% tra = full(data_vc2p.grad.balance.parcellation.tra(1:2,1:151));
% data_vc2p.trial = tra * cellrowselect(data_cmb.trial, 1:151);
% data_vc2p.trialinfo = data_cmb.trialinfo;
% data_vc2p.time = data_cmb.time;

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [0 100];
cfg.pad    = 1;
cfg.tapsmofrq = 5;
freq2 = ft_freqanalysis(cfg, data_vc2p);

cfg = [];
cfg.trials = find(freq2.trialinfo == 0);
bsl2 = ft_freqdescriptives(cfg, freq2);
cfg.trials = find(freq2.trialinfo == 1);
act2 = ft_freqdescriptives(cfg, freq2);

figure; plot(bsl2.freq, act2.powspctrm./bsl2.powspctrm - 1);

