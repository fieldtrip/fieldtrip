function test_data2bids

% WALLTIME 00:20:00
% MEM 5gb
% DEPENDENCY data2bids

%% Example with only behavioural data
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test1'))

cfg = [];
cfg.method                      = 'convert';
cfg.dataset                     = [];
cfg.presentationfile            = './sourcedata/A2002-2-MEG-MOUS-Aud.log';
cfg.writejson                   = 'replace';
cfg.bidsroot                    = fullfile(tempdir, 'bids');
cfg.sub                         = 'A2002';
cfg.ses                         = 'MEG';
cfg.datatype                    = 'events';
cfg.TaskName                    = 'language';
cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.dataset_description.Name    = 'test';
data2bids(cfg)


%% Example with a CTF dataset on disk that needs no conversion
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test2'))

cfg = [];
cfg.dataset                     = 'Subject01.ds';
cfg.TaskName                    = 'language';
cfg.meg.PowerLineFrequency      = 50;
cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
data2bids(cfg)

%% Example with a Neuroscan EEG dataset on disk that needs to be converted
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test3'))

cfg = [];
cfg.method                      = 'convert';
cfg.dataset                     = 'Subject1_MP.cnt';
cfg.datatype                    = 'eeg';
cfg.outputfile                  = 'sub-MP_task-visual_eeg.vhdr';
cfg.writejson                   = 'replace';
cfg.writetsv                    = 'replace';
cfg.TaskName                    = 'visual';
cfg.TaskDescription             = 'Visual response task';
cfg.Instructions                = 'Press as fast with your right index finger as you see the target appear on screen.';
cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
% give some extra information about two of the channels
cfg.channels.name = {'AF3', 'AF4'};
cfg.channels.status = {'bad', 'good'};
cfg.channels.impedance = [20e3 5e3];
data2bids(cfg)

%% Example with a BrainVision EEG dataset on disk that needs to be copied
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test4'))

cfg = [];
cfg.method                      = 'copy';
cfg.dataset                     = 'Subject1_MP.vhdr';
cfg.datatype                    = 'eeg';
cfg.bidsroot                    = tempname();
cfg.sub                         = '01';
cfg.ses                         = '02';
cfg.run                         = 3;
cfg.TaskName                    = 'visual';
cfg.TaskDescription             = 'Visual response task';
cfg.Instructions                = 'Press as fast with your right index finger as you see the target appear on screen.';
cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
% give some extra information about one of the channels
cfg.channels.name = 'AF3';
cfg.channels.status = 'bad';
cfg.channels.impedance = Inf;
data2bids(cfg)

%% Example with preprocessed EEG data in memory
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test5'))

cfg = [];
cfg.dataset                     = 'Subject1_MP.cnt';
cfg.bpfilter                    = 'yes';
cfg.bpfreq                      = [0.1 40];
data = ft_preprocessing(cfg);

cfg = [];
cfg.outputfile                  = 'sub-MP_task-visual_eeg.vhdr';
cfg.datatype                    = 'eeg'; % could also be ieeg or emg
cfg.writejson                   = 'yes';
cfg.writetsv                    = 'yes';
cfg.TaskName                    = 'visual';
cfg.TaskDescription             = 'Visual response task';
cfg.Instructions                = 'Press as fast with your right index finger as you see the target appear on screen.';
cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
data2bids(cfg, data)

%% Example with an anatomical MRI on disk that needs no conversion
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test6'))

cfg = [];
cfg.dataset                     = 'oostenveld_r.nii';
cfg.mri.dicomfile               = 'ERIVDBER_030731_R.OOSTERVELD.MR.PAUGAA_ANATOMICAL-3D.2.142.2003.7.31.11.19.16.140000.53832621.IMA';
cfg.mri.StationName             = 'Sonata';
cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.writejson                   = 'merge';
data2bids(cfg)

%% Example with realigned and resliced anatomical MRI data in memory
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test7'))

mri = ft_read_mri('dicom/ERIVDBER_030731_R.OOSTERVELD.MR.PAUGAA_ANATOMICAL-3D.2.56.2003.7.31.11.19.16.468000.53834351.IMA');

cfg = [];
cfg.outputfile                  = 'sub-RO_T1w.nii';
cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
data2bids(cfg, mri)

%% Example where the data is copied
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test8'));

cfg = [];
cfg.method                      = 'copy';
cfg.dataset                     = 'Eeglab_data.set';
cfg.outputfile                  = 'sub-EEG_task-attention_eeg.set';
cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
data2bids(cfg)

cfg = [];
cfg.method                      = 'copy';
cfg.dataset                     = 'jg_single_01raw.fif';
cfg.events                      = ft_read_event(cfg.dataset, 'chanindx', 123); % this has analog triggers, which require some extra attention
cfg.outputfile                  = 'sub-MEG_task-stimuluation_meg.fif';
cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
data2bids(cfg)

%% try some iEEG data with electrodes

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test9'));

% this loads data, elec_acpc_fr, and mri
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/human_ecog/SubjectUCI29/SubjectUCI29_data.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/human_ecog/SubjectUCI29/SubjectUCI29_elec_acpc_fr.mat'));
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/human_ecog/SubjectUCI29/SubjectUCI29_MR_acpc.nii'));

% make the data a bit smaller
cfg = [];
cfg.trials = 1;
data = ft_selectdata(cfg, data);

cfg = [];
cfg.resamplefs = 500;
data = ft_resampledata(cfg, data);

% add the electrode definition
data.elec = elec_acpc_fr;

cfg = [];
cfg.bidsroot = 'human_ecog';
cfg.sub = 'UCI29';
cfg.task = 'attention';
cfg.datatype = 'ieeg';
cfg.coordsystem.iEEGCoordinateSystem = 'ACPC';
cfg.coordsystem.iEEGCoordinateSystemDescription = 'electrodes were aligned with ACPC and projected on the cortex hull';
data2bids(cfg, data);

cfg = [];
cfg.bidsroot = 'human_ecog';
cfg.sub = 'UCI29';
cfg.datatype = 'T1w';
cfg.mri.StationName = 'bay3';
cfg.mri.MagneticFieldStrength = 3;
data2bids(cfg, mri);

%%
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test10'));

cfg = [];
cfg.dataset                     = 'sub-A2002/meg/sub-A2002_task-auditory_meg.ds';
cfg.writejson                   = 'replace';
cfg.presentationfile            = 'sourcedata/meg/A2002-2-MEG-MOUS-Aud.log';

cfg.writejson                   = 'replace';
cfg.writetsv                    = 'replace';

cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';

cfg.TaskName                    = 'Auditory language task';

cfg.trigger.eventtype           = 'UPPT001';
cfg.presentation.skip           = 'none';
cfg.trigger.eventvalue          = 3;
cfg.presentation.eventtype      = 'Nothing';
cfg.presentation.eventvalue     = '3 Audi*';

cfg.meg.DigitizedLandmarks      = true;
cfg.meg.DigitizedHeadPoints     = true;
cfg.meg.PowerLineFrequency      = 50;
cfg.meg.DewarPosition           = 'upright';
cfg.meg.SoftwareFilters         = 'n/a';

data2bids(cfg);

%%
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/test11'));

cfg = [];
cfg.dataset                     = 'sub-A2002/func/sub-A2002_task-auditory_bold.nii';
cfg.presentationfile            = 'sourcedata/mri_task/A2002-2-fMRI-MOUS-Aud.log';

cfg.writejson                   = 'merge';
cfg.writetsv                    = 'replace';

cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';

cfg.presentation.eventtype      = 'Pulse';
cfg.presentation.eventvalue     = [];
cfg.RepetitionTime              = 2;

data2bids(cfg);
