function test_beamformer_sam(dataset)
 
% MEM 1gb
% WALLTIME 00:20:00
% DEPENDENCY ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis ft_inverse_sam
% DATA public

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
cfg.sam.lambda          = '5%';
cfg.sam.projectmom      = 'no';
cfg.sam.keepfilter      = 'yes';
cfg.sam.keepori         = 'yes';
cfg.sam.projectnoise    = 'yes';

% default case, equivalent to the fixedori approach previously named 'moiseev'
cfg.sam.noisecov        = tlck_N.cov;
source_moiseev_ft       = ft_sourceanalysis(cfg, tlck_S);
ori_moiseev             = cat(2,source_moiseev_ft.avg.ori{:});

% if no noise covariance is specified, it is estimated using a multiple of the identity matrix.
% this is equivalent to the fixedori approach previously named 'gareth'
cfg.sam                 = rmfield(cfg.sam, 'noisecov');
source_gar_ft           = ft_sourceanalysis(cfg, tlck_S);
ori_gar                 = cat(2,source_gar_ft.avg.ori{:});

ori_All = [ori_moiseev; ori_gar];


%% estimates moms per trials

cfg                     = [];
cfg.method              = 'sam';
cfg.grid                = grid;
cfg.headmodel           = vol;
cfg.grad                = data_signal.grad;
cfg.sam.lambda          = '5%';
cfg.sam.keepfilter      = 'yes';
cfg.sam.keepori         = 'yes';
cfg.rawtrial            = 'yes';
cfg.sam.keepmom         = 'yes';  % saves the time series
cfg.sam.keeptrials      = 'yes'; %'no' or 'yes
cfg.sam.keepleadfield   = 'no';
cfg.sam.keepfilter      = 'no';

% default case, equivalent to the fixedori approach previously named 'moiseev'
cfg.sam.noisecov        = tlck_N.cov;
cfg.grid.filter         = source_moiseev_ft.avg.filter;
source_moiseev_mom      = ft_sourceanalysis(cfg, tlck_S);
 

% if no noise covariance is specified, it is estimated using a multiple of the identity matrix.
% this is equivalent to the fixedori approach previously named 'gareth'
cfg.sam                 = rmfield(cfg.sam, 'noisecov');
cfg.grid.filter         = source_gar_ft.avg.filter;
source_gar_mom          = ft_sourceanalysis(cfg, tlck_S);


% make ft data structure with original simulated signals
source_ori_mom = source_moiseev_mom;
for ith = 1:numel(source_ori_mom.trial)
  for iSrc = 1:numel(inx_in)
    source_ori_mom.trial(ith).mom{inx_in(iSrc)} = data_sinSig_c(iSrc,:);
  end
end

%% make FT dataset

source_moiseev_ft        = struct();
source_moiseev_ft.fsample= data_signal.fsample;
source_moiseev_ft.time   = {source_moiseev_mom.time};
source_moiseev_ft.label  = cellfun( @num2str, num2cell([1:numel(inx_in)]'), 'un', 0);
nsrc                     = numel(source_moiseev_ft.label);
source_gar_ft            = source_moiseev_ft;
source_ori_ft            = source_moiseev_ft;

for d = 1:ntrials
  
  source_moiseev_ft.trial{d} = cat(1,source_moiseev_mom.trial(d).mom{inx_in}); % data_ori(:,:,d);%  Noise_dat_o(:,:,d);%
  source_gar_ft.trial{d}     = cat(1,source_gar_mom.trial(d).mom{inx_in}); % data_ori(:,:,d);%  Noise_dat_o(:,:,d);%
  source_ori_ft.trial{d}     = cat(1,source_ori_mom.trial(d).mom{inx_in});
  
  source_moiseev_ft.time{d}  = source_moiseev_ft.time{1};
  source_gar_ft.time{d}      = source_gar_ft.time{1};
  source_ori_ft.time{d}      = source_ori_ft.time{1};
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
freq_moiseev     = ft_freqanalysis(cfg, source_moiseev_ft);
freq_gar         = ft_freqanalysis(cfg, source_gar_ft);
freq_ori         = ft_freqanalysis(cfg, source_ori_ft);

conn = {  'coh'  'plv' };
connout = { 'cohspctrm'  'plvspctrm' };

c = 1;
cfg           = [];
cfg.method    = conn{c};
      
cohm1          = ft_connectivityanalysis(cfg, freq_moiseev);
cohm2          = ft_connectivityanalysis(cfg, freq_gar);
cohm4          = ft_connectivityanalysis(cfg, freq_ori);
    
cfg = [];
cfg.parameter   = connout{1};
ft_connectivityplot(cfg, cohm4)
