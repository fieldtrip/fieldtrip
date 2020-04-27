function test_pull1377
% using addpath /home/common/matlab/fieldtrip
% MEM ?
% WALLTIME ?

% DEPENDENCY ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis 

% this function creates a set of input-structures to be used for testing source analysis steps prior to inverse solution

% struct of the script:
% 1. get the data (MEG and EEG)
% 2. preprocess the data (MEG and EEG): (timelock with covariance for lcmv beamforming and minimumnormestimate, 
%                                        timelock without keeptrials for dipolefitting,
%                                        freq for dics beamforming and pcc beamforming)
% 3. create leadfield: 3.a: internally computed leadfield (prepare volume conductor first and leadfield afterwards (MEG localspheres, EEG single sphere))
%                      3.b: load externally computed leadfield and convert to appropriate format
% 4. perform source analysis (MEG and EEG): LCMV beamforming, MNE, DICS, PCC, dipolefitting

%% 1. get the data
% for MEG data + sensor info
dataname = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat');
load(dataname);
datameg = data; 
clear data

% for EEG data + channel info
cfg=[];
cfg.dataset = '/home/common/matlab/fieldtrip/data/ftp/tutorial/preprocessing_erp/s04.eeg';
hdr        = ft_read_header(cfg.dataset);
event      = ft_read_event(cfg.dataset);
EVsample   = [event.sample]';
EVvalue    = {event.value}';
Word = find(strcmp('S141', EVvalue)==1);

begsample = EVsample(Word) - 10;
endsample = EVsample(Word) + 500;
offset = -100*ones(size(endsample));
cfg.trl = [begsample endsample offset];

dataeeg = ft_preprocessing(cfg);
% Fix EEG channel names
cfg=[];
cfg.layout = 'biosemi64.lay';
lay = ft_prepare_layout(cfg);
dataeeg.label = lay.label(1:end-2);
clear lay

%
% !! create worse-case scenario, whereby order and nr of chans don't match across inputs (grad/elec and leadfields)
% Step 1: remove 2-3 random chans from both MEG and EEG raw data
cfg=[];
cfg.channel = randperm(length(datameg.label)-3);
datameg = ft_selectdata(cfg, datameg);
cfg.channel = randperm(length(dataeeg.label)-2);
dataeeg = ft_selectdata(cfg, dataeeg);
% Setp 2: shuffle order of labels
%datameg.label = datameg.label(randperm(length(datameg.label)));
dataeeg.label = dataeeg.label(randperm(length(dataeeg.label)));

%% 2. preprocess the data
% create timelock structure with covariance for lcmv beamforming and minimumnormestimate, 
% and timeloch without keeptrials fro dipolefitting
% for MEG
cfg  = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
cfg.channel    = 'MEG';
MEG_tlck = ft_timelockanalysis(cfg, datameg);

cfg  = [];
cfg.covariance = 'yes';
% cfg.keeptrials = 'yes'; %if this is not commented, the .avg field necessary in dipolefitting is missing
cfg.channel    = 'MEG';
MEG_tlck_df = ft_timelockanalysis(cfg, datameg);


% for EEG 
cfg  = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
EEG_tlck = ft_timelockanalysis(cfg, dataeeg);

% for EEG (dimension problems in dipolefitting!!)
cfg  = [];
cfg.covariance = 'yes';
% cfg.keeptrials = 'yes';
EEG_tlck_df = ft_timelockanalysis(cfg, dataeeg);

% create freq structure for dics beamforming and pcc beamforming
% for MEG,
cfg  = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.tapsmofrq = 4;
cfg.foilim = [0 20];
cfg.channel = 'MEG';
MEG_freq = ft_freqanalysis(cfg, datameg);

% for EEG
cfg  = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.tapsmofrq = 4;
cfg.foilim = [0 20];
EEG_freq = ft_freqanalysis(cfg, dataeeg);

%% 3. create leadfield 
% 3.a internally

% get volume conductor model
% for MEG, localsphere
load('/home/common/matlab/fieldtrip/data/test/latest/vol/Subject01vol_localspheres.mat')
vol_localsphere     = vol;
clear vol

% for EEG, singlesphere
vol_singlesphere.o = [0 0 2];
vol_singlesphere.r = 8;
vol_singlesphere.c = 1;
vol_singlesphere.unit = 'cm';

% create leadfield (internally to ft)
% for MEG
cfg      = [];
cfg.grad = datameg.grad; 
cfg.headmodel = vol_localsphere;
cfg.channel = 'MEG';
cfg.resolution = 1.5;
gridmeg = ft_prepare_leadfield(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create 2D grid - is this necessary?? Not working
%[pnt, tri] = mesh_sphere(162);
%pnt   = pnt*(vol.orig.MEG_Sphere.RADIUS-1.5);
%shift = [vol.orig.MEG_Sphere.ORIGIN_X vol.orig.MEG_Sphere.ORIGIN_Y vol.orig.MEG_Sphere.ORIGIN_Z];
%pnt   = pnt+repmat(shift,[size(pnt,1) 1]);
%grid2.pnt = pnt;
%grid2.tri = tri;
%grid2.inside = 1:size(grid2.pnt,1);
%grid2.outside = [];

%cfg      = [];
%cfg.grad = datameg.grad;
%cfg.headmodel = vol;
%cfg.sourcemodel = grid2;
%cfg.channel = 'MEG';
%grid2 = ft_prepare_leadfield(cfg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for EEG
elecname = '/home/common/matlab/fieldtrip/template/electrode/standard_1020.elc';
elec = ft_read_sens(elecname, 'senstype', 'eeg');

cfg      = [];
cfg.grad = elec; % why not dataeeg.elec?
cfg.headmodel = vol_singlesphere;
cfg.channel = 'EEG';
cfg.resolution = 1.5;
grideeg = ft_prepare_leadfield(cfg); % inconsistent with dataeeg.labels

%% 3.b externally (mimic externally created leadfields)
% for MEG
leadf = randn(3000, length(datameg.grad)); % here we are assuming that the externally computed leadfield has the structure Ndip x Nsens
ext_leadfield_meg.unit  = 'cm';
ext_leadfield_meg.leadfielddimord = '{pos}_chan_ori';

% the leadfield
Ndip = size(leadf,1)/3;
ext_leadfield_meg.leadfield = cell(1,Ndip);

for d =1:Ndip
    ext_leadfield_meg.leadfield{d} = [leadf(d,:); leadf(d+Ndip,:); leadf(d+2*Ndip,:)]';
end

ext_leadfield_meg.inside = ones(size(ext_leadfield_meg.leadfield));
ext_leadfield_meg.pos = randn(length(ext_leadfield_meg.leadfield),3);


% for EEG
leadf = randn(3000, length(elec.label)); % here we are assuming that the externally computed leadfield has the structure Ndip x Nsens
ext_leadfield_eeg.unit  = 'cm';
ext_leadfield_eeg.leadfielddimord = '{pos}_chan_ori';

% the leadfield
Ndip = size(leadf,1)/3;
ext_leadfield_eeg.leadfield = cell(1,Ndip);

for d =1:Ndip
    ext_leadfield_eeg.leadfield{d} = [leadf(d,:); leadf(d+Ndip,:); leadf(d+2*Ndip,:)]';
end

ext_leadfield_eeg.inside = ones(size(ext_leadfield_eeg.leadfield));
ext_leadfield_eeg.pos = randn(length(ext_leadfield_eeg.leadfield),3);

%% 4. perform source analysis
% for MEG
% do LCMV beamforming
cfg            = [];
cfg.method     = 'lcmv';
cfg.lcmv.keepleadfield = 'yes';
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.keepcov       = 'yes';
cfg.lcmv.lambda        = '5%';
cfg.sourcemodel       = gridmeg;
cfg.headmodel  = vol_localsphere;
sourcelcmv3d1  = ft_sourceanalysis(cfg, MEG_tlck);
% project through computed filter
cfg.rawtrial    = 'yes';
cfg.sourcemodel        = gridmeg;
cfg.sourcemodel.filter = sourcelcmv3d1.avg.filter;
ft_sourceanalysis(cfg, MEG_tlck);


% do MNE 
cfg = [];
cfg.method   = 'mne';
cfg.mne.keepleadfield = 'yes';
cfg.mne.keepfilter = 'yes';
cfg.mne.lambda     = 1e4;
cfg.headmodel = vol_localsphere;
cfg.sourcemodel = gridmeg;
sourcemne3d1 = ft_sourceanalysis(cfg, MEG_tlck); 
cfg.rawtrial    = 'yes';
cfg.sourcemodel        = gridmeg;
cfg.sourcemodel.filter = sourcemne3d1.avg.filter;
ft_sourceanalysis(cfg, MEG_tlck);


% do DICS
cfg = [];
cfg.method = 'dics';
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.lambda        = '5%';
cfg.frequency = 10;
cfg.headmodel = vol_localsphere;
cfg.sourcemodel = gridmeg;
sourcedics3d1 = ft_sourceanalysis(cfg, MEG_freq);

cfg.rawtrial    = 'yes';
cfg.sourcemodel        = gridmeg;
cfg.sourcemodel.filter = sourcedics3d1.avg.filter;
ft_sourceanalysis(cfg, MEG_freq);

% do PCC 
cfg = [];
cfg.method = 'pcc';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.keepleadfield = 'yes';
cfg.pcc.keepcsd       = 'yes';
cfg.pcc.keepmom       = 'yes';
cfg.pcc.lambda        = '5%';
cfg.frequency = 10;
cfg.headmodel = vol_localsphere;
cfg.sourcemodel = gridmeg;
ft_sourceanalysis(cfg, MEG_freq);

% Only one inverse solution suffices since fixes to ft_sourceanalysis need to happen before inverse solution
% do DICS for EEG
cfg = [];
cfg.method = 'dics';
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.lambda        = '5%';
cfg.frequency = 10;
cfg.headmodel = vol_singlesphere;
cfg.sourcemodel = grideeg;
cfg.elec = elec;
sourcedics3d1 = ft_sourceanalysis(cfg, EEG_freq);
%% %%%%%
%% ERROR
cfg.rawtrial    = 'yes';
cfg.sourcemodel        = grideeg;
cfg.sourcemodel.filter = sourcedics3d1.avg.filter;
ft_sourceanalysis(cfg, EEG_freq);

% do dipolefit
% for MEG
cfg = [];
cfg.numdipoles    = 1;                              
cfg.headmodel     = vol_localsphere;                     
cfg.grid          = gridmeg;                     
cfg.nonlinear     = 'no';                           
cfg.grad          = MEG_tlck_df.grad;                        
cfg.latency       = 0.025;                         
ft_dipolefitting(cfg,MEG_tlck_df);

%% %%%%%%
%% ERROR
%% %%%%%%
% for EEG
cfg = [];
cfg.numdipoles    =  1;             
cfg.headmodel     = vol_singlesphere;         
cfg.grid          = grideeg;                 
cfg.elec          = elec;                         
cfg.latency       = 0.025;                          
ft_dipolefitting(cfg,EEG_tlck_df); %problems with dimensions!
% Is this because ft_dipolefitting doesn't behave when elec/leadfield/data have difference labels and order?

%% %%%%%%%%%%%%%%%%%%%%%%
% Previous ones work, what about externally generated leadfield:
%% %%%%%%%%%%%%%%%%%%%%%%
% Works if cfg.headmodel is provided
%% Used to work, no longer working
cfg = [];
cfg.method = 'dics';
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.lambda        = '5%';
cfg.frequency = 10;
cfg.headmodel = vol_localsphere;
cfg.sourcemodel = ext_leadfield_meg;
sourcedics3d1 = ft_sourceanalysis(cfg, MEG_freq);

% Does NOT work if cfg.headmodel is missing
cfg = [];
cfg.method = 'dics';
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.lambda        = '5%';
cfg.frequency = 10;
%cfg.headmodel = []; % vol_localsphere;
cfg.sourcemodel = ext_leadfield_meg;
sourcedics3d1 = ft_sourceanalysis(cfg, MEG_freq);







