function inspect_qsubcellfun

% MEM 3gb
% WALLTIME 00:31:04

% TEST inspect_qsubcellfun
% TEST qsubcellfun qsubfeval qsubget

if isempty(which('qsubcellfun'))
  [ftver, ftpath] = ft_version;
  addpath(fullfile(ftpath, 'qsub'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this test script is based on http://www.fieldtriptoolbox.org/tutorial/distributedcomputing

cd(dccnpath('/home/common/matlab/fieldtrip/data'));
timreq = 20*60;
memreq = 2*1024^3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subjectlist = {
  'Subject01.ds'
  'Subject02.ds'
  'Subject03.ds'
  'Subject04.ds'
  };
 
conditionlist = {
  'FC'
  'FIC'
  'IC'
  };
 
triggercode = [
  9
  3
  5
  ];
 
% start with a new and empty configuration
cfg = {};
 
for subj=1:4
  for cond=1:3
    cfg{subj,cond}                      = [];
    cfg{subj,cond}.dataset              = subjectlist{subj};
    cfg{subj,cond}.trialdef.prestim     = 1;
    cfg{subj,cond}.trialdef.poststim    = 2;
    cfg{subj,cond}.trialdef.eventtype   = 'backpanel trigger';
    cfg{subj,cond}.trialdef.eventvalue  = triggercode(cond);
  end
end
 
cfg = qsubcellfun(@ft_definetrial, cfg, 'timreq', timreq, 'memreq', memreq);
 
% this extends the previous configuration
for subj=1:4
  for cond=1:3
    cfg{subj,cond}.channel              = {'MEG', '-MLP31', '-MLO12'};
    cfg{subj,cond}.demean               = 'yes';
    cfg{subj,cond}.baselinewindow       = [-0.2 0];
    cfg{subj,cond}.lpfilter             = 'yes';
    cfg{subj,cond}.lpfreq               = 35;
  end
end
 
data = qsubcellfun(@ft_preprocessing, cfg, 'timreq', timreq, 'memreq', memreq);
 
% start with a new and empty configuration
cfg = {};
 
for subj=1:4
  for cond=1:3
    % timelockanalysis does not require any non-default settings
    cfg{subj,cond} = [];
  end
end
 
timelock = qsubcellfun(@ft_timelockanalysis, cfg, data, 'timreq', timreq, 'memreq', memreq);
 
% from here on we won't process the data in parallel any more
% average each condition over all subjects
cfg = [];
avgFC  = ft_timelockgrandaverage(cfg, timelock{:,1});
avgFIC = ft_timelockgrandaverage(cfg, timelock{:,2});
avgIC  = ft_timelockgrandaverage(cfg, timelock{:,3});
 
% cfg = [];
% cfg.layout = 'CTF151.lay';
% ft_multiplotER(cfg, avgFC, avgFIC, avgIC);
 
% cfg = [];
% cfg.channel = {'MLC33', 'MLC43', 'MLP11', 'MLP12', 'MLP13', 'MLP33', 'MLP34', 'MLT14', 'MLT15', 'MLT25'}
% ft_singleplotER(cfg, avgFC, avgFIC, avgIC);

