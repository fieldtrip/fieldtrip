function [Cf, Cr, Pr, Ntrials, cfg] = prepare_freq_matrices(cfg, freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that converts a freq structure into Cf, Cr and Pr
% this is used in FT_SOURCEANALYSIS
%
% This function returns data matrices with a channel order that is consistent
% with the original channel order in the data.
%
% The order of the channels in the output data is according to the input cfg.channel,
% which therefore must be specified as a cell-array with actual labels, not as an
% input like 'all' that still needs to be interpreted by FT_CHANNELSELECTION.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2015, Jan-Mathijs Schoffelen
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

% set the defaults
cfg = ft_checkconfig(cfg, 'deprecated', 'dicsfix');
cfg = ft_checkconfig(cfg, 'forbidden', {'frequency', 'latency'});
cfg = ft_checkconfig(cfg, 'required', 'channel'); % this should be a full list

cfg.keeptrials = ft_getopt(cfg, 'keeptrials', 'yes');
cfg.rawtrial   = ft_getopt(cfg, 'rawtrial', 'no');
cfg.refchan    = ft_getopt(cfg, 'refchan');

keeptrials = istrue(cfg.keeptrials) || istrue(cfg.rawtrial);

Cf = [];
Cr = [];
Pr = [];

if startsWith(freq.dimord, 'rpt')
  Ntrials = size(freq.cumtapcnt,1);
else
  Ntrials = 1;
end

% create a square csd-matrix, if necessary
hasfull = false;
if isfield(freq, 'crsspctrm')
  dimtok  = tokenize(getdimord(freq, 'crsspctrm'),'_');
  hasfull = sum(strcmp(dimtok, 'chan'))==2;
end
if ~hasfull
  if keeptrials
    freq = ft_checkdata(freq, 'cmbrepresentation', 'full');
  else
    freq = ft_checkdata(freq, 'cmbrepresentation', 'fullfast');
    Ntrials = 1;
  end
end

% extract the csd-matrix for the channels-of-interest, order them according to cfg.channel
[dum, chanindx] = match_str(cfg.channel, freq.label);

% update the cfg
cfg.channel = freq.label(chanindx);
if startsWith(freq.dimord, 'rpt')
  Cf = freq.crsspctrm(:,chanindx,chanindx,:,:);
else
  Cf = freq.crsspctrm(chanindx,chanindx,:,:);
end

if ~isempty(cfg.refchan)
  refindx = match_str(freq.label, cfg.refchan);
  if isempty(refindx)
    ft_error('the requested reference channel is not found in the data');
  end
  if startsWith(freq.dimord, 'rpt')
    Cr = freq.crsspctrm(:,chanindx,refindx,:,:);
    Pr = freq.crsspctrm(:,refindx,refindx,:,:);
  else
    Cr = freq.crsspctrm(chanindx,refindx,:,:);
    Pr = freq.crsspctrm(refindx,refindx,:,:);
  end
end  

% do a sanity check on the cross-spectral-density matrix
if any(isnan(Cf(:)))
  ft_error('The cross-spectral-density matrix is not complete');
end
if any(isnan(Cr(:)))
  ft_error('The cross-spectral-density with the reference channel is not complete');
end
