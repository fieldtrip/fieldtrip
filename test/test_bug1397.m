function test_bug1397

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preprocessing ft_appenddata

% the following code was obtained from http://www.fieldtriptoolbox.org/tutorial/coherence
% on Wed Mar 28 15:36:40 CEST 2012

% find the interesting epochs of data
cfg = [];
% MODIFICATION, use trialfun handle and other path to the data
cfg.trialfun                  = @trialfun_left;
cfg.dataset                   = dccnpath('/home/common/matlab/fieldtrip/data/SubjectCMC.ds');
cfg = ft_definetrial(cfg);

% MODIFICATION, use only 10 trials
cfg.trl = cfg.trl(1:10,:);

% MODIFICATION, the following should not affect the problem
%
% % detect EOG artifacts in the MEG data
% cfg.continuous                = 'yes';
% cfg.artfctdef.eog.padding     = 0;
% cfg.artfctdef.eog.bpfilter    = 'no';
% cfg.artfctdef.eog.detrend     = 'yes';
% cfg.artfctdef.eog.hilbert     = 'no';
% cfg.artfctdef.eog.rectify     = 'yes';
% cfg.artfctdef.eog.cutoff      = 2.5;
% cfg.artfctdef.eog.interactive = 'no';
% cfg = ft_artifact_eog(cfg);
%
% % detect jump artifacts in the MEG data
% cfg.artfctdef.jump.interactive = 'no';
% cfg.padding                    = 5;
% cfg = ft_artifact_jump(cfg);
%
% % detect muscle artifacts in the MEG data
% cfg.artfctdef.muscle.cutoff      = 8;
% cfg.artfctdef.muscle.interactive = 'no';
% cfg = ft_artifact_muscle(cfg);
%
% % reject the epochs that contain artifacts
% cfg.artfctdef.reject          = 'complete';
% cfg = ft_rejectartifact(cfg);

% preprocess the MEG data
cfg.demean                    = 'yes';
cfg.dftfilter                 = 'yes';
cfg.channel                   = {'MEG'};
cfg.continuous                = 'yes';
meg = ft_preprocessing(cfg);

cfg              = [];
cfg.dataset      = meg.cfg.dataset;
cfg.trl          = meg.cfg.trl;
cfg.continuous   = 'yes';
cfg.demean       = 'yes';
cfg.dftfilter    = 'yes';
cfg.channel      = {'EMGlft' 'EMGrgt'};
cfg.hpfilter     = 'yes';
cfg.hpfreq       = 10;
cfg.rectify      = 'yes';
emg = ft_preprocessing(cfg);

% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1397
% the reported problem in fieldtrip-20120302 was
% 
% ??? Error using ==> ft_appenddata at 266
% there is a difference in the time axes of the input data

data = ft_appenddata([], meg, emg);
assert(numel(data.label)==153);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to avoid external dependencies of this test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trl = trialfun_left(cfg)

% read in the triggers and create a trial-matrix
% consisting of 1-second data segments, in which
% left ECR-muscle is active.

event = ft_read_event(cfg.dataset);
trig  = [event(find(strcmp('backpanel trigger', {event.type}))).value];
indx  = [event(find(strcmp('backpanel trigger', {event.type}))).sample];

%left-condition
sel = [find(trig==1028):find(trig==1029)];

trig = trig(sel);
indx = indx(sel);

trl = [];
for j = 1:length(trig)-1
  trg1 = trig(j);
  trg2 = trig(j+1);
  if trg1<=100 && trg2==2080
    trlok      = [[indx(j)+1:1200:indx(j+1)-1200]' [indx(j)+1200:1200:indx(j+1)]'];
    trlok(:,3) = [0:-1200:-1200*(size(trlok,1)-1)]';
    trl        = [trl; trlok];
  end
end
