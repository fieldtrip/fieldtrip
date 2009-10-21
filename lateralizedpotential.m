function lrp = lateralizedpotential(cfg, avgL, avgR);

% LATERALIZEDPOTENTIAL computes lateralized potentials such as the LRP
%
% Use as
%   [lrp] = lateralizedpotential(cfg, avgL, avgR)
%
% where the input datasets should come from TIMELOCKANALYSIS
% and the configuration should contain
%   cfg.channelcmb = Nx2 cell array
%
% An example channelcombination containing the homologous channels 
% in the 10-20 standard system is
%    cfg.channelcmb = {'Fp1'   'Fp2'
%                      'F7'    'F8'
%                      'F3'    'F4'
%                      'T7'    'T8'
%                      'C3'    'C4'
%                      'P7'    'P8'
%                      'P3'    'P4'
%                      'O1'    'O2'}
%
% The lateralized potential is computed on combinations of channels and
% not on indivudual channels. However, if you want to make a topographic
% plot with e.g. MULTIPLOTER, you can replace the output lrp.label
% with lrp.plotlabel.
%
% The concept for the LRP was introduced approximately simultaneously in the
% following two papers
% - M. G. H. Coles. Modern mind-brain reading - psychophysiology,
%   physiology, and cognition. Psychophysiology, 26(3):251-269, 1988.
% - R. de Jong, M. Wierda, G. Mulder, and L. J. Mulder. Use of
%   partial stimulus information in response processing. J Exp Psychol
%   Hum Percept Perform, 14:682-692, 1988.
% and it is discussed in detail on a technical level in
% - R. Oostenveld, D.F. Stegeman, P. Praamstra and A. van Oosterom.
%   Brain symmetry and topographic analysis of lateralized event-related
%   potentials. Clin Neurophysiol. 114(7):1194-202, 2003.
%
% See also LATERALIZEDFIELD, TIMELOCKANALYSIS, MULTIPLOTER

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: lateralizedpotential.m,v $
% Revision 1.6  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.5  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.4  2005/06/24 08:26:43  roboos
% restructured the help, changed some whitespace
%
% Revision 1.3  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.2  2005/01/25 08:13:40  roboos
% fixed small error with brackets
%
% Revision 1.1  2004/12/16 11:19:05  roboos
% new implementation
%

fieldtripdefs

if ~isfield(cfg, 'channelcmb'), 
  cfg.channelcmb = {
    'Fp1'   'Fp2'
    'F7'    'F8'
    'F3'    'F4'
    'T7'    'T8'
    'C3'    'C4'
    'P7'    'P8'
    'P3'    'P4'
    'O1'    'O2'
    };
end

if ~all(avgL.time==avgR.time)
  error('timeaxes are not the same');
end

% start with an empty output structure
lrp.label     = {};
lrp.plotlabel = {};
lrp.avg       = [];
lrp.time      = avgL.time;

% compute the lateralized potentials
Nchan = size(cfg.channelcmb);
for i=1:Nchan
  C3R = strmatch(cfg.channelcmb{i,1}, avgR.label);
  C4R = strmatch(cfg.channelcmb{i,2}, avgR.label);
  C3L = strmatch(cfg.channelcmb{i,1}, avgL.label);
  C4L = strmatch(cfg.channelcmb{i,2}, avgL.label);
  if ~isempty(C3R) && ~isempty(C4R) && ~isempty(C3L) && ~isempty(C4L)
    lrp.label{end+1}     = sprintf('%s/%s', cfg.channelcmb{i,1}, cfg.channelcmb{i,2});
    lrp.plotlabel{end+1} = cfg.channelcmb{i,1};
    erpC3L = avgL.avg(C3L,:);
    erpC4L = avgL.avg(C4L,:);
    erpC3R = avgR.avg(C3R,:);
    erpC4R = avgR.avg(C4R,:);
    lrp.avg(end+1,:) = 1/2 * ((erpC3R - erpC4R) + (erpC4L - erpC3L));
  end
end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: lateralizedpotential.m,v 1.6 2008/09/22 20:17:43 roboos Exp $';
% remember the configuration details of the input data
cfg.previous = [];
try, cfg.previous{1} = avgL.cfg; end
try, cfg.previous{2} = avgR.cfg; end
% remember the exact configuration details in the output 
lrp.cfg = cfg;

