function [cfg] = definetrial(cfg);

% DEFINETRIAL defines the trials, i.e. the pieces of data that will be read
% in for preprocessing. Trials are defined by their begin and end sample
% in the data file and each trial has an offset that defines where the
% relative t=0 point (usually the point of the trigger) is for that trial.
%
% Use as
%   [cfg] = definetrial(cfg)
% where the configuration structure should contain either
%   cfg.trialdef   = structure with details of trial definition, see below
%   cfg.trialfun   = function name, see below
%
% A call to DEFINETRIAL results in the trial definition "trl" being added
% to the output configuration structure. The trials are defined according
% to the triggers, trials or other events in the data, or from a
% user-specified Matlab function which returns "trl".
%
% The trial definition "trl" is an Nx3 matrix, N is the number of trials.
% The first column contains the sample-indices of the begin of each trial 
% relative to the begin of the raw data, the second column contains the 
% sample-indices of the end of each trial, and the third column contains 
% the offset of the trigger with respect to the trial. An offset of 0 
% means that the first sample of the trial corresponds to the trigger. A 
% positive offset indicates that the first sample is later than the trigger, 
% a negative offset indicates that the trial begins before the trigger.
%
% Simple trial definitions (e.g. based on a trigger alone) are supported by
% DEFINETRIAL itself. For this, the general and data format independent way
% of handling trials is by relying on the READ_FCDC_EVENT function to
% collect all event information (such as triggers) from your dataset and
% select trials based on those events. This is implemented in DEFINETRIAL as
%   cfg.trialdef.eventtype  = 'string'
%   cfg.trialdef.eventvalue = number, string or list with numbers or strings
%   cfg.trialdef.prestim    = number, latency in seconds (optional)
%   cfg.trialdef.poststim   = number, latency in seconds (optional)
%
% If you specify cfg.trialdef.eventtype  = '?' a list with the events in your 
% data file will be displayed on screen.
%
% However, there are also many other complex ways in which you can define
% data pieces of interest, for example based on a conditional sequence of
% events (e.g. stimulus trigger followed by a correct response). For those
% cases, a general mechanism has been implemented through which you can
% supply your own trial-defining function, the 'trialfun'.
%
% This 'trialfun' is a string containing the name of a function that you
% have to write yourself. The function should take the cfg-structure as
% input and should give a Nx3 matrix in the same format as "trl" as the
% output. You can add extra custom fields to the configuration structure to
% pass as arguments to your own trialfun. Furthermore, inside the trialfun
% you can use the READ_EVENT function to get the event information
% from your data file.
%
% See also PREPROCESSING, READ_HEADER, READ_DATA, READ_EVENT

% Undocumented local options:
% cfg.datafile
% cfg.dataset
% cfg.event
% cfg.trl
% cfg.version

% Copyright (c) 2003, Robert Oostenveld, F.C. Donders Centre
%
% $Log: definetrial.m,v $
% Revision 1.57  2009/10/07 12:42:16  roevdmei
% changed reference to older FCDC read functions
%
% Revision 1.56  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.55  2009/01/14 11:29:55  sashae
% temporarily disabled previous revision
%
% Revision 1.54  2009/01/13 10:14:58  sashae
% changed handling of the output cfg: now the cfg also has cfg.previous fields,
% similar to data.cfg.previous. this way the output of definetrial and the
% artifact functions is kept separately from subsequent preprocessing steps
%
% Revision 1.53  2008/10/10 14:41:22  sashae
% replaced call to dataset2files with checkconfig
%
% Revision 1.52  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.51  2007/07/27 12:36:24  roboos
% updated a comment in the code
%
% Revision 1.50  2007/05/06 08:54:11  roboos
% do not add version info to cfg, will be done by preprocessing
%
% Revision 1.49  2006/09/18 07:41:50  roboos
% in case of eventtype=? give a gentle information message (instead of error)
%
% Revision 1.48  2006/06/20 16:25:58  ingnie
% updated documentation
%
% Revision 1.47  2006/05/30 20:57:15  roboos
% updated documentation
%
% Revision 1.46  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.45  2006/03/29 10:34:02  ingnie
% updated documentation
%
% Revision 1.44  2005/11/04 16:06:15  ingnie
% fixed some typos in help
%
% Revision 1.43  2005/10/04 16:12:25  roboos
% changed exception handling of nargout for matlab65 into try-catch
%
% Revision 1.42  2005/09/07 10:18:31  roboos
% in matlab65 (R13) it does not find the functions in the private directory, whereas in matlab70 it does
% hence, made the specific evaluation of the nargout function dependent on the matlab version using "which"
%
% Revision 1.41  2005/09/02 14:24:14  roboos
% rigorously cleaned up the code
% moved all dataformat specific code into separate subfunctions
% added defaults for cfg.trialfun if not specified, depending on datatype and settings
% removed old CVS log history
%
% Revision 1.40  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.39  2005/02/16 15:18:25  roboos
% cleaned up the documentation, removed cfg.trialfile as option (was not supported)
% cleaned up the handling of events, retain event if already present in cfg (just like trl)
% removed unclear warning about old v.s. new style of configuration
% improved fprintf feedback on number of events and trials
%

fieldtripdefs

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'dataset2files', {'yes'});

if ~isfield(cfg, 'trl') && (~isfield(cfg, 'trialfun') || isempty(cfg.trialfun))
  % try to determine a standard trial function for the various data formats
  % most of these trial functions are old-fashioned and do not use the
  % event information that is returned by READ_FCDC_EVENT

  % determine the trials from the general events structure
  if isfield(cfg, 'trialdef') && isfield(cfg.trialdef, 'eventtype')
    cfg.trialfun = 'trialfun_general';
    % determine the trials from the events in the EEProbe file
  elseif isfield(cfg, 'trialdef') && isfield(cfg.trialdef, 'trgfile') && filetype(cfg.trialdef.trgfile, 'eep_trg')
    cfg.trialfun = 'trialfun_eeprobe_cnt';
    % determine the trials from the average EEProbe file
  elseif filetype(cfg.datafile, 'eep_avr')
    cfg.trialfun = 'trialfun_eeprobe_avr';
    % determine the trials from the events in the BrainVision file
  elseif isfield(cfg, 'trialdef') && isfield(cfg.trialdef, 'trgfile') && filetype(cfg.trialdef.trgfile, 'brainvision_vmrk')
    cfg.trialfun = 'trialfun_brainvision';
    % determine the trials from the Neuromag file
  elseif filetype(cfg.datafile, 'neuromag_fif')
    cfg.trialfun = 'trialfun_neuromag';
    % determine the trials from the triggers or class-definitions in the epoched CTF file
  elseif filetype(cfg.dataset, 'ctf_ds') && isfield(cfg, 'trialdef') &&  (isfield(cfg.trialdef, 'includeTrigger') || isfield(cfg.trialdef, 'excludeTrigger') || isfield(cfg.trialdef, 'includeConditions') || isfield(cfg.trialdef, 'excludeConditions'))
    cfg.trialfun = 'trialfun_ctf_epoched';
    % determine the trials from the triggers in the continuous CTF file
  elseif filetype(cfg.dataset, 'ctf_ds') && isfield(cfg, 'trialdef') && ~(isfield(cfg.trialdef, 'includeTrigger') || isfield(cfg.trialdef, 'excludeTrigger') || isfield(cfg.trialdef, 'includeConditions') || isfield(cfg.trialdef, 'excludeConditions') || isfield(cfg.trialdef, 'eventtype'))
    cfg.trialfun = 'trialfun_ctf_continuous';
    % determine the trials from the epoched NeuroScan file
  elseif filetype(cfg.datafile, 'ns_eeg')
    cfg.trialfun = 'trialfun_neuroscan_eeg';
    % determine the trials from the events in the continuous NeuroScan file
  elseif isfield(cfg, 'trialdef') && filetype(cfg.datafile, 'ns_cnt')
    cfg.trialfun = 'trialfun_neuroscan_cnt';
  else
    error('no trialfunction specified, see DEFINETRIAL for help');
  end
end

% create the trial definition for this dataset and condition
if isfield(cfg, 'trl')
  % the trial definition is already part of the configuration
  fprintf('retaining exist trial definition\n');
  trl = cfg.trl;
  if isfield(cfg, 'event')
    fprintf('retaining exist event information\n');
    event = cfg.event;
  else
    event = [];
  end
elseif isfield(cfg, 'trialfun')
  % evaluate the user-defined function that gives back the trial definition
  fprintf('evaluating trialfunction ''%s''\n', cfg.trialfun);
  % determine the number of outpout arguments of the user-supplied trial function
  try
    % the nargout function in Matlab 6.5 and older does not work on function handles
    num = nargout(cfg.trialfun);
  catch
    num = 1;
  end
  if num==1
    % the user-defined function only gives back the trial definition
    trl   = feval(cfg.trialfun, cfg);
    event = [];
  else
    % the user-defined function also gives back detailed information about
    % conditions, reaction time or any other information
    [trl, event] = feval(cfg.trialfun, cfg);
  end
else
  error('no trialfunction specified, see DEFINETRIAL for help');
end

if isfield(cfg, 'trialdef') && isfield(cfg.trialdef, 'eventtype') && strcmp(cfg.trialdef.eventtype, '?')
  % give a gentle message instead of an error
  fprintf('no trials have been defined yet, see DEFINETRIAL for further help\n');
elseif size(trl,1)<1
  error('no trials were defined, see DEFINETRIAL for help');
end

% add the new trials and events to the output configuration
fprintf('found %d events\n', length(event));
cfg.event = event;
fprintf('created %d trials\n', size(trl,1));
cfg.trl = trl;

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add information about the version of this function to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i1] = dbstack;
  cfg.version.name = st(i1);
end
cfg.version.id = '$Id: definetrial.m,v 1.57 2009/10/07 12:42:16 roevdmei Exp $';

% % remember the exact configuration details in the output
% cfgtmp = cfg;
% cfg = [];
% try cfg.trl        = cfgtmp.trl;        end
% try cfg.dataset    = cfgtmp.dataset;    end
% try cfg.datafile   = cfgtmp.datafile;   end
% try cfg.headerfile = cfgtmp.headerfile; end
% cfg.previous = cfgtmp;
