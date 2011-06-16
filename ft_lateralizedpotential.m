function lrp = ft_lateralizedpotential(cfg, avgL, avgR);

% FT_LATERALIZEDPOTENTIAL computes lateralized potentials such as the LRP
%
% Use as
%   [lrp] = ft_lateralizedpotential(cfg, avgL, avgR)
%
% where the input datasets should come from FT_TIMELOCKANALYSIS
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
% plot with e.g. FT_MULTIPLOTER, you can replace the output lrp.label
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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_TIMELOCKANALYSIS, FT_MULTIPLOTER

% Copyright (C) 2004, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();

% set the defaults
if ~isfield(cfg, 'inputfile'),  cfg.inputfile                   = [];    end
if ~isfield(cfg, 'outputfile'), cfg.outputfile                  = [];    end

hasdata = nargin>1;
if ~isempty(cfg.inputfile) % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    avgL=loadvar(cfg.inputfile{1}, 'data'); % read first element as avgL from array inputfile
    avgR=loadvar(cfg.inputfile{2}, 'data'); % read second element as avgR from array inputfile
  end
end

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
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();

% remember the configuration details of the input data
cfg.previous = [];
try, cfg.previous{1} = avgL.cfg; end
try, cfg.previous{2} = avgR.cfg; end

% remember the exact configuration details in the output 
lrp.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', lrp); % use the variable name "data" in the output file
end

