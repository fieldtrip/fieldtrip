function test_bug3359

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivityanalysis ft_topoplotER 

% first create some data
%--------------------------------------------------------
% make 3 channels with no direct link between 1 and 2
cfg             = [];
cfg.ntrials     = 200;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';
cfg.params(:,:,1) = [ 0.8 0   0; 
                      0   0.9 0.5;
                      0.4 0   0.5];
cfg.params(:,:,2) = [-0.5    0  0; 
                        0 -0.8  0; 
                        0    0 -0.2];
cfg.noisecov      = [0.3 0 0;
                       0 1 0;
                       0 0 0.2];

data_sim          = ft_connectivitysimulation(cfg);

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'));

cfg         = [];
cfg.channel = 'MEG';
data        = ft_selectdata(cfg, data);

mixing_matrix = blkdiag(ones(71,1), ones(80,1), 1);
data.trial = cell(1,numel(data_sim.trial));
for k = 1:numel(data_sim.trial)
  data.trial{k} = mixing_matrix*data_sim.trial{k} + randn(152,200)./2;
end
data.time           = data_sim.time;
data.label{end+1,1} = 'refchannel';
  
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [0 80];
cfg.taper  = 'hanning';
cfg.pad    = 1;
freq = ft_freqanalysis(cfg, data);
freq = ft_checkdata(freq, 'cmbstyle', 'fullfast');

freq_sim = ft_freqanalysis(cfg, data_sim);
freq_sim = ft_checkdata(freq_sim, 'cmbstyle', 'fullfast');

cfg = [];
cfg.method = 'psi';
cfg.bandwidth = 5;
psi = ft_connectivityanalysis(cfg, freq);

% compute this for reference
psi_sim = ft_connectivityanalysis(cfg, freq_sim);

cfg = [];
cfg.parameter = 'psispctrm';
figure;ft_connectivityplot(cfg, psi_sim);

cfg = [];
cfg.parameter = 'psispctrm';
cfg.layout    = 'CTF151_helmet.mat';
cfg.refchannel = 'MLC11';
cfg.directionality = 'inflow';
ft_topoplotER(cfg, psi);
cfg.directionality = 'outflow';
ft_topoplotER(cfg, psi);

cfg.refchannel = 'refchannel';
ft_topoplotER(cfg, psi);

% linearly indexed channel combinations, three flavours
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'powandcsd';
cfg.foilim = [0 80];
cfg.taper  = 'hanning';
cfg.pad    = 1;
cfg.channelcmb = ft_channelcombination({'MEG' 'refchannel'}, data.label, 0, 0);
freq0 = ft_freqanalysis(cfg, data);
cfg.channelcmb = ft_channelcombination({'MEG' 'refchannel'}, data.label, 0, 1);
freq1 = ft_freqanalysis(cfg, data);
cfg.channelcmb = ft_channelcombination({'MEG' 'refchannel'}, data.label, 0, 2);
freq2 = ft_freqanalysis(cfg, data);

cfg           = [];
cfg.method    = 'psi';
cfg.bandwidth = 5;
psi0          = ft_connectivityanalysis(cfg, freq0);
psi1          = ft_connectivityanalysis(cfg, freq1);
psi2          = ft_connectivityanalysis(cfg, freq2);

cfg           = [];
cfg.parameter = 'psispctrm';
cfg.layout    = 'CTF151_helmet.mat';
cfg.refchannel = 'refchannel';
ft_topoplotER(cfg, psi0);
ft_topoplotER(cfg, psi1);
try
  ft_topoplotER(cfg, psi2);
catch
  fprintf('running ft_topoplotER with ''chancmb'' without cfg.directionality leads to a problem with psi2\n');
end
cfg.directionality = 'inflow';
ft_topoplotER(cfg, psi0);
try
  ft_topoplotER(cfg, psi1);
catch
  fprintf('error psi1 plotting with ''inflow''\n');
end
ft_topoplotER(cfg, psi2);

cfg.directionality = 'outflow';
try
  ft_topoplotER(cfg, psi0);
catch
  fprintf('error psi0 plotting with ''outflow''\n');
end
ft_topoplotER(cfg, psi1);
ft_topoplotER(cfg, psi2);
