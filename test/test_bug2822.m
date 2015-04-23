function test_bug2822

% TEST test_bug2822
% TEST ft_sourceanalysis minimumnormestimate 


% implement frequency domain MNE reconstruction


%% Load data (dataFIC)
path_to_load = ('/home/common/matlab/fieldtrip');
load(dccnpath([path_to_load, '/data/test/dataFIC.mat']))

%% Frequency analysis

cfg              = [];
cfg.output       = 'powandcsd';
cfg.channel      = 'MEG';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 18;                          
freq = ft_freqanalysis(cfg, dataFIC);

%% Fake Headmodel and leadfield
load(dccnpath([path_to_load, '/template/headmodel/standard_singleshell.mat']))

% use 'icosahedron' private function to generate the mash
mesh = [];
[mesh.pnt, mesh.tri] = icosahedron642();
mesh.pnt = 5*mesh.pnt - repmat([ 0 3 -1.5],size(mesh.pnt,1),1) ;

%% Leadfield

cfg = [];
cfg.grad = dataFIC.grad;                      % sensor positions
cfg.channel = {'MEG', '-MLP31', '-MLO12'};   % the used channels
cfg.grid.pos = mesh.pnt;              % source points
cfg.grid.inside = 1:size(mesh.pnt,1); % all source points are inside of the brain
cfg.vol = vol;                               % volume conduction model
leadfield = ft_prepare_leadfield(cfg);

%% Source analysis MNE

cfg = [];
cfg.method = 'mne';
cfg.frequency = 18;
cfg.grid   = leadfield;
cfg.vol    = vol;
cfg.lambda = 0.001;
source_freq_mne = ft_sourceanalysis(cfg, freq);

%% Source analysis RV

cfg = [];
cfg.method = 'rv';
cfg.frequency = 18;
cfg.grid   = leadfield;
cfg.vol    = vol;
cfg.lambda = 0.001;
source_freq_rv = ft_sourceanalysis(cfg, freq);

%% Source analysis MUSIC

cfg = [];
cfg.method = 'music';
cfg.frequency = 18;
cfg.grid   = leadfield;
cfg.vol    = vol;
cfg.music.numcomponent = 1;
cfg.lambda = 0.001;
source_freq_music = ft_sourceanalysis(cfg, freq);




