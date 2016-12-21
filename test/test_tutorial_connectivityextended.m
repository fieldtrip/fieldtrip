function test_tutorial_connectivityextended

% WALLTIME 00:45:00
% MEM 3gb

% TEST test_tutorial_connectivity
% TEST ft_connectivityanalysis ft_connectivitysimulation ft_freqanalysis ft_connectivityplot ft_mvaranalysis

% This is the first section of the connectivity tutorial, which
% starts with an MVAR model and then uses parametric and nonparametric
% spectral decomposition for coherence and granger

% See also test_tutorial_connectivity2 and test_tutorial_connectivity3

global ft_default;
ft_default.feedback = 'no';

%% simulate data
cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';
cfg.params(:,:,1) = [ 0.8    0    0 ;
  0  0.9  0.5 ;
  0.4    0  0.5];

cfg.params(:,:,2) = [-0.5    0    0 ;
  0 -0.8    0 ;
  0    0 -0.2];

cfg.noisecov      = [ 0.3    0    0 ;
  0    1    0 ;
  0    0  0.2];
data              = ft_connectivitysimulation(cfg);

figure
plot(data.time{1}, data.trial{1})
legend(data.label)
xlabel('time (s)')

cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
ft_databrowser(cfg, data);

%% mvaranalysis
cfg         = [];
cfg.order   = 5;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data);

%% freqanalysis 1
cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);

%% freqanalysis 2
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, data);

%% connectivityanalysis
cfg           = [];
cfg.method    = 'coh';
coh           = ft_connectivityanalysis(cfg, freq);
cohm          = ft_connectivityanalysis(cfg, mfreq);

%% visualisation
cfg           = [];
cfg.parameter = 'cohspctrm';
ft_connectivityplot(cfg, coh, cohm);

%% do the same for granger
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, mfreq);

cfg           = [];
cfg.parameter = 'grangerspctrm';
ft_connectivityplot(cfg, granger);

figure
for row=1:3
  for col=1:3
    subplot(3,3,(row-1)*3+col);
    plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)))
    ylim([0 1])
  end
end

%% do the virtual channel stuff
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/source_coh_lft.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/source_diff.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/data_cmb.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/sourcemodel.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/hdm.mat'));

[maxval, maxcohindx] = max(source_coh_lft.avg.coh);
source_coh_lft.pos(maxcohindx, :)
[maxval, maxpowindx] = max(source_diff.avg.pow);
source_diff.pos(maxpowindx, :)

cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = 'MEG';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
tlock                 = ft_timelockanalysis(cfg, data_cmb);

cfg              = [];
cfg.method       = 'lcmv';
cfg.headmodel    = hdm;
cfg.grid.pos     = sourcemodel.pos([maxcohindx maxpowindx], :);
cfg.grid.inside  = true(2,1);
cfg.grid.unit    = sourcemodel.unit;
cfg.lcmv.keepfilter = 'yes';
source_idx       = ft_sourceanalysis(cfg, tlock);

beamformer_lft_coh = source_idx.avg.filter{1};
beamformer_gam_pow = source_idx.avg.filter{2};

chansel = ft_channelselection('MEG', data_cmb.label); % find MEG sensor names
chansel = match_str(data_cmb.label, chansel);         % find MEG sensor indices

coh_lft_data = [];
coh_lft_data.label = {'coh_lft_x', 'coh_lft_y', 'coh_lft_z'};
coh_lft_data.time = data_cmb.time;
gam_pow_data = [];
gam_pow_data.label = {'gam_pow_x', 'gam_pow_y', 'gam_pow_z'};
gam_pow_data.time = data_cmb.time;
for i=1:length(data_cmb.trial)
  coh_lft_data.trial{i} = beamformer_lft_coh * data_cmb.trial{i}(chansel,:);
  gam_pow_data.trial{i} = beamformer_gam_pow * data_cmb.trial{i}(chansel,:);
end

cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
%ft_databrowser(cfg, gam_pow_data);

visualTimeseries = cat(2, gam_pow_data.trial{:});
motorTimeseries = cat(2, coh_lft_data.trial{:});
[u1, s1, v1] = svd(visualTimeseries, 'econ');
[u2, s2, v2] = svd(motorTimeseries, 'econ');

virtualchanneldata = [];
virtualchanneldata.label = {'visual', 'motor'};
virtualchanneldata.time = data_cmb.time;

for k = 1:length(data_cmb.trial)
  virtualchanneldata.trial{k}(1,:) = u1(:,1)' * beamformer_gam_pow * data_cmb.trial{k}(chansel,:);
  virtualchanneldata.trial{k}(2,:) = u2(:,1)' * beamformer_lft_coh * data_cmb.trial{k}(chansel,:);
end

% select the two EMG channels
cfg = [];
cfg.channel = 'EMG';
emgdata = ft_selectdata(cfg, data_cmb);

% combine the virtual channel with the two EMG channels
cfg = [];
combineddata = ft_appenddata(cfg, virtualchanneldata, emgdata);

%% compute the spectral decomposition
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel    = {'visual' 'motor' 'EMGlft' 'EMGrgt'};
freq    = ft_freqanalysis(cfg, combineddata);

cfg = [];
cfg.method = 'coh';
coherence = ft_connectivityanalysis(cfg, freq);

cfg = [];
cfg.zlim = [0 0.25];
figure
ft_connectivityplot(cfg, coherence);

