function test_data2bids

% WALLTIME 00:20:00
% MEM 5gb

% DEPENDENCY data2bids

%%
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/mous'));

cfg = [];
cfg.dataset                     = 'sub-A2002/meg/sub-A2002_task-auditory_meg.ds';
cfg.meg.writesidecar            = 'replace';
cfg.presentationfile            = 'sourcedata/meg/A2002-2-MEG-MOUS-Aud.log';

cfg.meg.writesidecar            = 'replace';
cfg.channels.writesidecar       = 'replace';
cfg.events.writesidecar         = 'replace';
cfg.coordsystem.writesidecar    = 'replace';

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
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/mous'));

cfg = [];
cfg.dataset                     = 'sub-A2002/func/sub-A2002_task-auditory_bold.nii';
cfg.presentationfile            = 'sourcedata/mri_task/A2002-2-fMRI-MOUS-Aud.log';

cfg.mri.writesidecar            = 'merge';
cfg.events.writesidecar         = 'replace';

cfg.InstitutionName             = 'Radboud University';
cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
cfg.InstitutionAddress          = 'Kapittelweg 29, 6525 EN, Nijmegen, The Netherlands';

cfg.presentation.eventtype      = 'Pulse';
cfg.presentation.eventvalue     = [];
cfg.RepetitionTime              = 2;

data2bids(cfg);

%%
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/data2bids/eeglab'));

% the following uses a dataset that originates from the EEGLAB website
% wget --no-check-certificate https://sccn.ucsd.edu/mediawiki/images/9/9c/Eeglab_data.set

age = [11  96  nan 77  82  87  18 40  26  80];
sex = {'f' [] 'f' 'f' 'f' 'm' 'm' 'm' 'm' 'm'};

sub = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'};
ses = {'pre', 'post'};
run = {1, 2};

for subindx=1:numel(sub)
  for sesindx=1:numel(ses)
    for runindx=1:numel(run)
      
      cfg = [];
      cfg.dataset   = 'sourcedata/Eeglab_data.set';
      cfg.datatype  = 'eeg';
      
      cfg.participants.age = age(subindx);
      cfg.participants.sex = sex{subindx};
      
      cfg.scans.acq_time = datestr(now, 'yyyy-mm-ddThh:MM:SS'); % RFC3339
      
      cfg.InstitutionName             = 'Radboud University';
      cfg.InstitutionalDepartmentName = 'Donders Institute for Brain, Cognition and Behaviour';
      
      cfg.bidsroot  = 'bids';
      cfg.sub       = sub{subindx};
      cfg.ses       = ses{sesindx};
      cfg.run       = run{runindx};
      data2bids(cfg);
      
    end % for run
  end % for ses
end % for sub


