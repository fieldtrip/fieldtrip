function test_bug1142

% MEM 1500mb
% WALLTIME 00:10:00

% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1142
%
% TEST ft_definetrial ft_preprocessing ft_rejectvisual ft_rejectartifact
% TEST ft_artifact_clip ft_artifact_eog ft_artifact_muscle ft_artifact_zvalue ft_artifact_ecg ft_artifact_jump ft_artifact_threshold

% the script is not yet finished, so at this moment it should not proceed to give an error
return

datadir = dccnpath('/home/common/matlab/fieldtrip/data');

artifact_type = {
  'clip'
  'eog'
  'muscle'
  'zvalue'
  'ecg'
  'jump'
  'threshold'
  };

order    = perms(1:length(artifact_type));
numorder = size(order,1);

% there are too many permutations to test them all, make a random subselection
order    = order(randperm(numorder),:);
numorder = 15;

% hopefully the defaults will do for most artifact detection functions
artfctdef_clip = [];

artfctdef_eog = [];

artfctdef_muscle = [];

artfctdef_zvalue = [];
artfctdef_zvalue.channel = 'MEG';
artfctdef_zvalue.cutoff = 4;
artfctdef_zvalue.interactive = 'no';

artfctdef_ecg = [];

artfctdef_jump = [];

artfctdef_threshold = [];
artfctdef_threshold.bpfreq = [0.5 40];

% collect the artifact configurations in the same order
artifact_cfg = {
  artfctdef_clip
  artfctdef_eog
  artfctdef_muscle
  artfctdef_zvalue
  artfctdef_ecg
  artfctdef_jump
  artfctdef_threshold
  };

cfg = [];
cfg.dataset = fullfile(datadir, 'ArtifactMEG.ds');
cfg.trialdef.triallength = 1;   % duration in seconds
cfg.trialdef.ntrials     = 30;  % number of trials
cfg = ft_definetrial(cfg);

cfg.continuous = 'yes';
cfg.checkconfig = 'loose';

data = ft_preprocessing(cfg);

cfg.trl = cfg.trl(11:20,:);  % select a subset for artifact detection, stay away from the edge of the data

cfg_ini = cfg;

for i=1:numorder
  sequence = order(i,:);
  
  cfg_v1 = cfg_ini;
  cfg_v2 = cfg_ini;
  
  for j=1:length(sequence)
    cfg_v1.artfctdef.(artifact_type{sequence(j)}) = artifact_cfg{sequence(j)};
    cfg_v1 = feval(['ft_artifact_' artifact_type{sequence(j)}], cfg_v1);
    
    cfg_v2.artfctdef.(artifact_type{sequence(j)}) = artifact_cfg{sequence(j)};
    cfg_v2 = feval(['ft_artifact_' artifact_type{sequence(j)}], cfg_v2, data);
  end
  
  cfg_all{i,1} = cfg_v1;
  cfg_all{i,2} = cfg_v2;
  
  cfg_v1 = ft_rejectartifact(cfg_v1, data);
  
  data_all{i,1} = ft_preprocessing(cfg_v1);
  data_all{i,2} = ft_rejectartifact(cfg_v2, data);
  
end

