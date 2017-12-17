function [lrp] = ft_lateralizedpotential(cfg, avgL, avgR)

% FT_LATERALIZEDPOTENTIAL computes lateralized potentials such as the
% lateralized readiness potential (LRP)
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
% To facilitate data-handling and distributed computing you can use
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
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar avgL avgR
ft_preamble provenance avgL avgR
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
avgL = ft_checkdata(avgL, 'datatype', 'timelock');
avgR = ft_checkdata(avgR, 'datatype', 'timelock');

% set the defaults
if ~isfield(cfg, 'channelcmb')
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

if ~isequal(avgL.time, avgR.time)
  ft_error('the time axes are not the same');
end

% start with an empty output structure
lrp.label     = {};
lrp.plotlabel = {};
lrp.avg       = [];
lrp.time      = avgL.time;

% add timelock signature
if isfield(avgL, 'dimord') && isfield(avgR, 'dimord')
    if ~strcmp(avgL.dimord, avgR.dimord)
        ft_error('The input data are of different dimord types');
    else
        lrp.dimord = avgL.dimord;
    end
else
    ft_error('''dimord'' not found. The function expects timelock data');
end

% compute the lateralized potentials
Nchan = size(cfg.channelcmb);
for i=1:Nchan
  % here the channel names "C3" and "C4" are used to clarify the
  % computation of the lateralized potential on all channel pairs
  C3R = strcmp(cfg.channelcmb{i,1}, avgR.label);
  C4R = strcmp(cfg.channelcmb{i,2}, avgR.label);
  C3L = strcmp(cfg.channelcmb{i,1}, avgL.label);
  C4L = strcmp(cfg.channelcmb{i,2}, avgL.label);
  if any(C3R) && any(C4R) && any(C3L) && any(C4L)
    lrp.label{end+1}     = sprintf('%s/%s', cfg.channelcmb{i,1}, cfg.channelcmb{i,2});
    lrp.plotlabel{end+1} = cfg.channelcmb{i,1};
    erpC3L = avgL.avg(C3L,:);
    erpC4L = avgL.avg(C4L,:);
    erpC3R = avgR.avg(C3R,:);
    erpC4R = avgR.avg(C4R,:);
    lrp.avg(end+1,:) = 1/2 * ((erpC3R - erpC4R) + (erpC4L - erpC3L));
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   avgL avgR
ft_postamble provenance lrp
ft_postamble history    lrp
ft_postamble savevar    lrp
