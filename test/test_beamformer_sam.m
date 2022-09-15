function test_beamformer_sam(dataset)
 
% MEM 3gb
% WALLTIME 00:20:00
% DEPENDENCY ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis ft_inverse_sam

% this function creates a set of source-structures to be used for testing

if nargin==0
  dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
end

% get volume conductor model
volname = fullfile(dataset, 'default.hdm');
vol     = ft_read_headmodel(volname);

% get data + sensor info
cfg         = [];
cfg.dataset = dataset;
cfg.channel = 'meg';
data        = ft_preprocessing(cfg);

%% get only the intertrial period

cfg = [];
cfg.latency = [-1 0];
data = ft_selectdata(cfg, data);

% create 3D grid
cfg      = [];
cfg.grad = data.grad;
cfg.headmodel = vol;
cfg.channel = 'MEG';
cfg.sourcemodel.resolution = 1.5;
cfg.reducerank = 3; % no rank reduction
grid = ft_prepare_leadfield(cfg);

inx_in = [1421, 1503]; %% source locations to be used

grid.inside = zeros(size(grid.inside));
grid.inside(inx_in) = 1;


%% make simulation signals
ntrials = numel(data.trial);
winlen = data.fsample;
dt = 1/winlen;
nsecs = abs(sum(data.time{1}([1,end])));
t = 0:dt:(nsecs);
f1 = 10;
f2 = 35;

x1 = sin(2*pi*f1*t) +sin(2*pi*f2*t)  +randn(size(t));
x2 = sin(2*pi*f1*t) +sin(2*pi*f2*t)  +randn(size(t));

data_sinSig_c = [x1;x2];
nsmaples = size(data_sinSig_c,2);
data_sinSig = repmat(data_sinSig_c,1,ntrials);

U3D = [1 0 0];
LF = cat(3, grid.leadfield{inx_in});
LFu =[];
for ith = 1:numel(inx_in)
    LFu(ith,:) =  LF(:,:,ith) *U3D';
end
nsensors = size(LFu,2);
ts_s = data_sinSig' * LFu;
ts_s = reshape(ts_s, nsmaples, ntrials, nsensors);
ts_s = permute(ts_s, [3 1 2]);

%% make sensor datasets

data_noise = data;
data_signal = data;

ts_n = cat(3, data_noise.trial{:});

% scale signal wrt noise
snr_ratio = sqrt(mean(ts_n(:).^2))/sqrt(mean(ts_s(:).^2));
ts_s = ts_s * snr_ratio;

SNR = 4;
ts_s = ts_s/SNR;

for ith = 1:numel(data_signal.trial)
  data_signal.trial{ith} =   ts_s(:,:,ith) + ts_n(:,:,ith);
end


% create timelock structure with covariance for lcmv and mne
cfg             = [];
cfg.covariance  = 'yes';
cfg.keeptrials  = 'yes';
cfg.channel     = 'MEG';
tlck_S          = ft_timelockanalysis(cfg, data_signal);
cfg.keeptrials  = 'no';
tlck_N          = ft_timelockanalysis(cfg, data_noise);

% do SAM beamforming
cfg            = [];
cfg.method              = 'sam';
cfg.grid                = grid;
cfg.headmodel           = vol;
cfg.grad                = data_signal.grad;
%cfg.sam.toi            = [.1 .2];
cfg.sam.lambda          = '5%';
cfg.sam.projectmom      = 'no';
cfg.sam.keepfilter      = 'yes';
cfg.sam.keepori         = 'yes';
cfg.sam.projectnoise    = 'yes';
cfg.sam.noisecov        = tlck_N.cov;
cfg.sam.fixedori        = 'moiseev';
source_new_ft           = ft_sourceanalysis(cfg, tlck_S);
ori_newImpl             = cat(2,source_new_ft.avg.ori{:});

cfg.sam.fixedori            = 'gareth';
source_gar_ft           = ft_sourceanalysis(cfg, tlck_S);
ori_gar                 = cat(2,source_gar_ft.avg.ori{:});

cfg.sam.fixedori            = 'robert';
source_rob_ft           = ft_sourceanalysis(cfg, tlck_S);
ori_rob                 = cat(2,source_rob_ft.avg.ori{:});

ori_All = [ori_newImpl; ori_gar; ori_rob];


%% estimates moms per trials

cfg                     = [];
cfg.method              = 'sam';
cfg.grid                = grid;
cfg.headmodel           = vol;
cfg.grad                = data_signal.grad;
%cfg.sam.toi            = [.1 .2];
cfg.sam.lambda          = '5%';
cfg.sam.keepfilter      = 'yes';
cfg.sam.keepori         = 'yes';
cfg.rawtrial            = 'yes';
cfg.sam.keepmom         = 'yes';  % saves the time series
cfg.sam.keeptrials      = 'yes'; %'no' or 'yes
cfg.sam.keepleadfield   = 'no';
cfg.sam.keepfilter      = 'no';
cfg.sam.fixedori        = 'moiseev';
cfg.sam.noisecov        = tlck_N.cov;
cfg.grid.filter         = source_new_ft.avg.filter;
source_new_mom          = ft_sourceanalysis(cfg, tlck_S);
 
cfg.sam.fixedori        = 'gareth';
cfg.grid.filter         = source_gar_ft.avg.filter;
source_gar_mom          = ft_sourceanalysis(cfg, tlck_S);

cfg.sam.fixedori        = 'robert';
source_rob_mom          = ft_sourceanalysis(cfg, tlck_S);

% make ft data structure with original simulated signals
source_ori_mom = source_new_mom;
for ith = 1:numel(source_ori_mom.trial)
  for iSrc = 1:numel(inx_in)
    source_ori_mom.trial(ith).mom{inx_in(iSrc)} = data_sinSig_c(iSrc,:);
  end
end

%% make FT dataset

source_new_ft            = struct();
source_new_ft.fsample    = data_signal.fsample;
source_new_ft.time       = {source_new_mom.time};
source_new_ft.label      = cellfun( @num2str, num2cell([1:numel(inx_in)]'), 'un', 0);
nsrc                     = numel(source_new_ft.label);
source_gar_ft            = source_new_ft;
source_rob_ft            = source_new_ft;
source_ori_ft            = source_new_ft;

for d = 1:ntrials
  
  source_new_ft.trial{d} = cat(1,source_new_mom.trial(d).mom{inx_in}); % data_ori(:,:,d);%  Noise_dat_o(:,:,d);%
  source_gar_ft.trial{d} = cat(1,source_gar_mom.trial(d).mom{inx_in}); % data_ori(:,:,d);%  Noise_dat_o(:,:,d);%
  source_rob_ft.trial{d} = cat(1,source_rob_mom.trial(d).mom{inx_in});
  source_ori_ft.trial{d} = cat(1,source_ori_mom.trial(d).mom{inx_in});
  
  source_new_ft.time{d} = source_new_ft.time{1};
  source_gar_ft.time{d} = source_gar_ft.time{1};
  source_rob_ft.time{d} = source_rob_ft.time{1};
  source_ori_ft.time{d} = source_ori_ft.time{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.foilim    = [1 50];
cfg.pad       =  'nextpow2';
cfg.tapsmofrq = 2;
freq_new         = ft_freqanalysis(cfg, source_new_ft);
freq_gar         = ft_freqanalysis(cfg, source_gar_ft);
freq_rob         = ft_freqanalysis(cfg, source_rob_ft);
freq_ori         = ft_freqanalysis(cfg, source_ori_ft);

conn = {  'coh'  'plv' };
connout = { 'cohspctrm'  'plvspctrm' };

c = 1;
cfg           = [];
cfg.method    = conn{c};
      
cohm1          = ft_connectivityanalysis(cfg, freq_new);
cohm2          = ft_connectivityanalysis(cfg, freq_gar);
cohm3          = ft_connectivityanalysis(cfg, freq_rob);
cohm4          = ft_connectivityanalysis(cfg, freq_ori);
    
cfg = [];
cfg.parameter   = connout{1};
ft_connectivityplot(cfg, cohm4)
