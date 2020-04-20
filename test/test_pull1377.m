function test_pull1377
% MEM 2gb
% WALLTIME 00:10:00

% DEPENDENCY ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis 

% this function creates a set of input-structures to be used for testing source analysis steps prior to inverse solution

%% get volume conductor model

% for MEG, localsphere
load('/home/common/matlab/fieldtrip/data/test/latest/vol/Subject01vol_localspheres.mat')
vol_localsphere     = vol;
clear vol

% for EEG, singlesphere?????
%volname = dccnpath('/home/common/matlab/fieldtrip/template/headmodel/standard_seg.mat'); %% is there an already segmented?
%cfg = [];
%cfg.method='singlesphere';
%vol_singlesphere = ft_prepare_headmodel(cfg, segmentedmri);

%%
%vol_singleshell = load('/home/common/matlab/fieldtrip/template/headmodel/standard_singleshell.mat');
vol_bem = load('/home/common/matlab/fieldtrip/template/headmodel/standard_bem.mat');

%% get MEG data + sensor info
dataname = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat');
load(dataname);
datameg = data; 
clear data

% get EEG data + channel info
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


% !! create worse-case scenario, whereby order and nr of chans don't match across inputs
% remove 2-3 random chans from both MEG and EEG raw data
cfg=[];
cfg.channel = randperm(length(datameg.label)-3);
datameg = ft_selectdata(cfg, datameg);
cfg.channel = randperm(length(dataeeg.label)-2);
dataeeg = ft_selectdata(cfg, dataeeg);


%% create 3D grid
% for MEG
cfg      = [];
cfg.grad = datameg.grad; %% check that this has remained the same size despite sensors removed
cfg.headmodel = vol_localsphere;
cfg.channel = 'MEG';
cfg.sourcemodel.resolution = 1.5;
gridmeg = ft_prepare_leadfield(cfg);

%% create 2D grid - is this necessary?? Not working
[pnt, tri] = mesh_sphere(162);
pnt   = pnt*(vol.orig.MEG_Sphere.RADIUS-1.5);
shift = [vol.orig.MEG_Sphere.ORIGIN_X vol.orig.MEG_Sphere.ORIGIN_Y vol.orig.MEG_Sphere.ORIGIN_Z];
pnt   = pnt+repmat(shift,[size(pnt,1) 1]);
grid2.pnt = pnt;
grid2.tri = tri;
grid2.inside = 1:size(grid2.pnt,1);
grid2.outside = [];

cfg      = [];
cfg.grad = datameg.grad;
cfg.headmodel = vol;
cfg.sourcemodel = grid2;
cfg.channel = 'MEG';
grid2 = ft_prepare_leadfield(cfg);
%%

% for EEG
elecname = '/home/common/matlab/fieldtrip/template/electrode/standard_1020.elc';
elec = ft_read_sens(elecname, 'senstype', 'eeg');

cfg      = [];
cfg.grad = elec;
cfg.headmodel = vol_singlesphere;
cfg.channel = 'EEG';
cfg.sourcemodel.resolution = 1.5;
grideeg = ft_prepare_leadfield(cfg);


%% mimic externally created leadfields
% for MEG
leadf = randn(442953, 151);
ext_leadfield.unit  = 'cm';
ext_leadfield.leadfielddimord = '{pos}_chan_ori';

% the leadfield
Ndip = size(leadf,1)/3;
ext_leadfield.leadfield = cell(1,Ndip);

for d =1:Ndip
    ext_leadfield.leadfield{d} = [leadf(d,:); leadf(d+Ndip,:); leadf(d+2*Ndip,:)]';
end

ext_leadfield.inside = ones(size(ext_leadfield.leadfield));
ext_leadfield.pos = randn(length(ext_leadfield.leadfield),3);


% for EEG - not sure what this would be... likely same structure, right?

%% sensor-level data
% create timelock structure with covariance for lcmv beamforming and
% minimumnormestimate
% for MEG
cfg  = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
cfg.channel    = 'MEG';
MEG_tlck = ft_timelockanalysis(cfg, datameg);

% for EEG 
cfg  = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
EEG_tlck = ft_timelockanalysis(cfg, dataeeg);

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


outputdir = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/source/meg'));

%% %%%%%%%%%%%%%%%%%%%%%
% create worse-case scenario, whereby order and nr of chans don't match across inputs
%% %%%%%%%%%%%%%%%%%%%%%
% sensor-level data has 2-3 less
% now remove other labels from vol, and change the order of labels




%% inverse solutions

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
%cfg.sourcemodel       = grid2;
%cfg.outputfile = 'ctf151_lcmv2d_avg';
%sourcelcmv2d1  = ft_sourceanalysis(cfg, MEG_tlck);

cfg.rawtrial    = 'yes';
cfg.sourcemodel        = gridmeg;
%cfg.outputfile  = fullfile(outputdir, 'ctf151_lcmv3d_trial');
cfg.sourcemodel.filter = sourcelcmv3d1.avg.filter;
ft_sourceanalysis(cfg, MEG_tlck);
%cfg.sourcemodel        = grid2;
%cfg.outputfile  = fullfile(outputdir, 'ctf151_lcmv2d_trial');
%cfg.sourcemodel.filter = sourcelcmv2d1.avg.filter;
%ft_sourceanalysis(cfg, MEG_tlck);

% do MNE 
cfg = [];
cfg.method   = 'mne';
cfg.mne.keepleadfield = 'yes';
cfg.mne.keepfilter = 'yes';
cfg.mne.lambda     = 1e4;
cfg.headmodel = vol_localsphere;
cfg.sourcemodel = gridmeg;
sourcemne3d1 = ft_sourceanalysis(cfg, tlck);
%cfg.sourcemodel = grid2;
%cfg.outputfile = fullfile(outputdir, 'ctf151_mne2d');
%sourcemne2d1 = ft_sourceanalysis(cfg, MEG_tlck);

%cfg.rawtrial    = 'yes';
%cfg.sourcemodel        = grid;
%cfg.sourcemodel.filter = sourcemne3d1.avg.filter;
%cfg.outputfile  = fullfile(outputdir, 'ctf151_mne3d_trial');
%ft_sourceanalysis(cfg, tlck);
%cfg.sourcemodel        = grid2;
%cfg.sourcemodel.filter = sourcemne2d1.avg.filter;
%cfg.outputfile  = fullfile(outputdir, 'ctf151_mne2d_trial');
%ft_sourceanalysis(cfg, MEG_tlck);

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
%cfg.outputfile = fullfile(outputdir, 'ctf151_dics3d_avg');
sourcedics3d1 = ft_sourceanalysis(cfg, MEG_freq);
%cfg.sourcemodel = grid2;
%cfg.outputfile = fullfile(outputdir, 'ctf151_dics2d_avg');
%sourcedics2d1 = ft_sourceanalysis(cfg, MEG_freq);

cfg.rawtrial    = 'yes';
cfg.sourcemodel        = gridmeg;
cfg.sourcemodel.filter = sourcedics3d1.avg.filter;
%cfg.outputfile  = fullfile(outputdir, 'ctf151_dics3d_trial');
ft_sourceanalysis(cfg, MEG_freq);
%cfg.sourcemodel        = grid2;
%cfg.sourcemodel.filter = sourcedics2d1.avg.filter;
%cfg.outputfile  = fullfile(outputdir, 'ctf151_dics2d_trial');
%ft_sourceanalysis(cfg, MEG_freq);


% do PCC - NOT TESTED
cfg = [];
cfg.method = 'pcc';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.keepleadfield = 'yes';
cfg.pcc.keepcsd       = 'yes';
cfg.pcc.keepmom       = 'yes';
cfg.pcc.lambda        = '5%';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.sourcemodel = grid;
cfg.outputfile = fullfile(outputdir, 'ctf151_pcc3d');
ft_sourceanalysis(cfg, freq);
cfg.sourcemodel = grid2;
cfg.outputfile = fullfile(outputdir, 'ctf151_pcc2d');
ft_sourceanalysis(cfg, MEG_freq);

% do dipolefit - NOT RUNNING
cfg = [];
cfg.numdipoles    = 1;                              %number of expected sources
cfg.headmodel     = vol_localsphere;                     %the head model
cfg.grid          = gridmeg;                      %the (precomputed) leadfield
cfg.nonlinear     = 'no';                           %only dipole scan
cfg.grad          = MEG_tlck.grad;                           %the sensor model
cfg.latency       = 0.025;                          %the latency of interest
%dipfit_meg    = ft_dipolefitting(cfg,MEG_avg);
dipfit_meg = ft_dipolefitting(cfg,MEG_tlck); % what should this input be?


% for EEG
cfg = [];
cfg.numdipoles    =  1;                             
cfg.headmodel     = headmodel;                  
cfg.grid          = leadfield;    
cfg.nonlinear     = 'no';                          
cfg.elec          = elec;                           
cfg.latency       = 0.025;                         
ft_dipolefitting(cfg,EEG_tlck);

