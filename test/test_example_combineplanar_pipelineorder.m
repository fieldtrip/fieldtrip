function test_example_combineplanar_pipelineorder

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_dipolesimulation ft_combineplanar ft_singleplotER ft_topoplotER

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

grad275 = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/megrealign/ctf275.mat'));

vol = [];
vol.r = 12;
vol.o = [0 0 4];

% in this example the dipole is carefully positioned on the same depth as
% where the dipole layer for the interpolation will be located. The center
% of the head is at [0 0 4], the radius is 12 cm, and the dipole layer that is
% specified in ft_megrealign is 2.5 cm shifted inward from the head surface

fsample = 500;
ntrial = 20;
smp = fsample*2; % 2 seconds epoch duration
time = -1:(1/fsample):1-(1/fsample); 
f = 10; % signal frequency content

cfg=[];
for k=1:ntrial;
  amp = randn() + 1; % vary signal amplitude trial-by-trial
  endsig = randi([35,50],1); % vary randomly the begining of the alpha desynchronization
  t = linspace(endsig-100,endsig,smp);
  env = (-1./(1+exp(-t/(rand()+5)))) + amp; % inverted sigmoid function to simulate alpha desynchronization
  
  phi = randn()*0.5; % vary signal phase trial-by-trial
  signal = cos(f*time*2*pi + phi);
  cfg.dip.signal{k} = (env .* signal) + randn(1,smp); %power sigmoid envelope * oscillation + gaussian random noise
end

cfg.dip.pos = [-1 -1 7];%dipole located on posterior sensors
cfg.fsample = fsample;
cfg.vol     = vol;
cfg.grad    = grad275;
data_axial  = ft_dipolesimulation(cfg);

%% synthetic planar computation
cfg              = [];
cfg.feedback     = 'no';
cfg.method       = 'template';
cfg.planarmethod = 'sincos';
cfg.channel      = {'MEG'};
cfg.trials       = 'all';
cfg.neighbours   = ft_prepare_neighbours(cfg, data_axial);
data_planar      = ft_megplanar(cfg,data_axial);

%% time-frequency analysis on axial and synthetic planar data
cfg            = [];
cfg.output     = 'pow';
cfg.channel    = 'MEG';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.foi        = 2:1:40;
cfg.toi        = [0:0.05:2];
cfg.trials     = 'all';
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
cfg.keeptrials = 'yes';
freq_planar    = ft_freqanalysis(cfg,data_planar);
freq_axial     = ft_freqanalysis(cfg,data_axial);

%% average the axial representation to compare later
freq_axial_avg = ft_freqdescriptives([],freq_axial);


%% Now let's test the influence of pipeline order on the espected value
%% A) combine planar on single trials and then average
cfg               = [];
cfg.combinemethod = 'sum';
freq_combined     = ft_combineplanar(cfg,freq_planar);

freq_combined_avg1 = ft_freqdescriptives([],freq_combined);

%% B) average single trials and then combine planar
freq_planar_avg = ft_freqdescriptives([],freq_planar);

cfg                = [];
cfg.combinemethod  = 'sum';
freq_combined_avg2 = ft_combineplanar(cfg,freq_planar_avg);


%%
cfg         = [];
cfg.layout  = 'CTF275.lay';
cfg.channel = 'MEG';
cfg.xlim    = [0.5 1];
cfg.ylim    = [8 12];

figure; 
subplot(321);ft_topoplotTFR(cfg, freq_axial_avg);title('axial');colorbar;
subplot(323);ft_topoplotTFR(cfg, freq_combined_avg1);title('combined-then-average');colorbar;
subplot(325);ft_topoplotTFR(cfg, freq_combined_avg2);title('average-then-combined');colorbar;

cfg         = [];
cfg.channel = {'MRO14', 'MRO24', 'MRO34', 'MRT27', 'MRT37'};%right axial cluster
subplot(322);ft_singleplotTFR(cfg, freq_axial);title('axial');colorbar;

cfg.channel = {'MRO11', 'MRO12', 'MRO14', 'MRO22', 'MRO24', 'MRO31', 'MRO32', 'MRO34', 'MRT27', 'MRT37'};%combined-planar cluster
subplot(324);ft_singleplotTFR(cfg, freq_combined_avg1);title('combined-then-average');colorbar;
subplot(326);ft_singleplotTFR(cfg, freq_combined_avg2);title('average-then-combined');colorbar;


%% Now let's going to play the same game but using Event-Related Fields (ERF)
%% A) combine planar on single trials and then average
cfg               = [];
cfg.combinemethod = 'sum';
data_combined     = ft_combineplanar(cfg,data_planar);

data_combined_avg1 = ft_timelockanalysis([],data_combined);


%% B) average single trials and then combine planar
data_axial_avg = ft_timelockanalysis([],data_axial);

cfg              = [];
cfg.feedback     = 'no';
cfg.method       = 'template';
cfg.planarmethod = 'sincos';
cfg.channel      = {'MEG'};
cfg.trials       = 'all';
cfg.neighbours   = ft_prepare_neighbours(cfg, data_axial_avg);
data_planar_avg  = ft_megplanar(cfg,data_axial_avg);

cfg                = [];
cfg.combinemethod  = 'sum';
data_combined_avg2 = ft_combineplanar(cfg,data_planar_avg);


%% C) average planar single trials and then combine planar
data_planar_avg = ft_timelockanalysis([],data_planar);

cfg                = [];
cfg.combinemethod  = 'sum';
data_combined_avg3 = ft_combineplanar(cfg,data_planar_avg);

%%
cfg         = [];
cfg.layout  = 'CTF275.lay';
cfg.channel = 'MEG';
cfg.xlim    = [0.5 1];

figure; 
subplot(421);ft_topoplotER(cfg, data_axial_avg);    title('axial');colorbar;
subplot(423);ft_topoplotER(cfg, data_combined_avg1);title('combined then average');colorbar;
subplot(425);ft_topoplotER(cfg, data_combined_avg2);title('average axial then planar+combined');colorbar;
subplot(427);ft_topoplotER(cfg, data_combined_avg3);title('average planar then combined');colorbar;

cfg = [];
cfg.channel = {'MRO14', 'MRO24', 'MRO34', 'MRT27', 'MRT37'};%right axial cluster
subplot(422);ft_singleplotER(cfg, data_axial_avg);title('axial');

cfg.channel = {'MRO11', 'MRO12', 'MRO14', 'MRO22', 'MRO24', 'MRO31', 'MRO32', 'MRO34', 'MRT27', 'MRT37'};%combined-planar cluster
subplot(424);ft_singleplotER(cfg, data_combined_avg1);title('combined then average');
subplot(426);ft_singleplotER(cfg, data_combined_avg2);title('average axial then planar+combined');
subplot(428);ft_singleplotER(cfg, data_combined_avg3);title('average planar then combined');

