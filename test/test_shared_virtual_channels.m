function failed_shared_virtual_channels

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_timelockanalysis ft_sourceanalysis ft_channelselection ft_databrowser

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/hdm.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/source_diff.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/source_coh_lft.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer_extended/data_cmb.mat'));

[maxval, maxcohindx] = max(source_coh_lft.avg.coh);
source_coh_lft.pos(maxcohindx, :)

assert(isalmostequal(source_coh_lft.pos(maxcohindx, :), [2.8000   -0.8000   7.2000], 'reltol', 0.001), 'coherence peak location not what it used to be!');

[maxval, maxpowindx] = max(source_diff.avg.pow);
source_diff.pos(maxpowindx, :)

assert(isalmostequal(source_diff.pos(maxpowindx, :), [0.4000   -8.8000    2.4000], 'reltol', 0.001), 'gamma power peak location not what it used to be!');

cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = 'MEG';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
tlock                 = ft_timelockanalysis(cfg, data_cmb);


cfg             = [];
cfg.method      = 'lcmv';
cfg.headmodel   = hdm;
cfg.sourcemodel.pos    = source_diff.pos([maxcohindx maxpowindx], :);
cfg.sourcemodel.inside = 1:size(cfg.sourcemodel.pos, 1);
cfg.sourcemodel.outside = [];
cfg.lcmv.keepfilter  = 'yes';
source_idx      = ft_sourceanalysis(cfg, tlock);

%% this is the old way of doing it.
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

%% this is the new way of doing it
cfg = [];
cfg.pos = source_idx.pos;
cfg.method = 'none';
data_vc = ft_virtualchannel(cfg, data_cmb, source_idx);
assert(isalmostequal(data_vc.trial{1}(4:6,:),gam_pow_data.trial{1}(1:3,:), 'reltol', 1e-9));

cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
ft_databrowser(cfg, gam_pow_data);
