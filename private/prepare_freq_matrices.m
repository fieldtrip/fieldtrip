function [Cf, Cr, Pr, Ntrials, cfg] = prepare_freq_matrices(cfg, freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that converts a freq structure into Cf, Cr and Pr
% this is used in sourecanalysis
%
% This function returns data matrices with a channel order that is consistent
% with the original channel order in the data.
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
if ~isfield(cfg, 'keeptrials'), cfg.keeptrials = 1;     end
if ~isfield(cfg, 'refchan'),    cfg.refchan    = [];    end
if ~isfield(cfg, 'rawtrial'),   cfg.rawtrial   = [];    end

keeptrials = istrue(cfg.keeptrials) || istrue(cfg.rawtrial);

Cf = [];
Cr = [];
Pr = [];

tok = tokenize(freq.dimord, '_');
if any(strcmp(tok, 'rpttap'))
  Ntrials = size(freq.cumtapcnt,1);
elseif any(strcmp(tok, 'rpt'))
  Ntrials = size(freq.cumtapcnt,1);
else
  Ntrials = 1;
end

% select from the frequency dimension
if any(strcmp(tok, 'freq')),
  % select the frequency of interest
  tmpcfg             = [];
  tmpcfg.frequency   = cfg.frequency;
  tmpcfg.avgoverfreq = 'yes';
  freq               = ft_selectdata(tmpcfg, freq);

  % update the cfg
  cfg.frequency      = freq.freq;
end

% select from the time dimension
if any(strcmp(tok, 'time')),
  % select the latency of interest for time-frequency data
  tmpcfg         = [];
  tmpcfg.latency = cfg.latency;
  tmpcfg.avgovertime = 'yes';
  freq           = ft_selectdata(tmpcfg, freq);
  
  % update the cfg
  cfg.latency    = freq.time;
end  

% create a square csd-matrix, if necessary
hasfull = false;
if isfield(freq, 'crsspctrm')
	dimtok  = tokenize(getdimord(freq, 'crsspctrm'),'_');
	hasfull = sum(strcmp(dimtok, 'chan'))==2;
end
if ~hasfull,
	if keeptrials,
		freq = ft_checkdata(freq, 'cmbrepresentation', 'full');
	else
		freq = ft_checkdata(freq, 'cmbrepresentation', 'fullfast');
		Ntrials = 1;
	end
end
tok = tokenize(freq.dimord, '_');

% extract the csd-matrix for the channels-of-interest
[dum, chanindx] = match_str(cfg.channel, freq.label);

% update the cfg
cfg.channel     = freq.label(chanindx);
if any(strncmp(tok, 'rpt', 3)),
  Cf = freq.crsspctrm(:,chanindx,chanindx);
else
  Cf = freq.crsspctrm(chanindx,chanindx);
end

if isfield(cfg, 'refchan') && ~isempty(cfg.refchan)
  refindx = match_str(freq.label, cfg.refchan);
  if isempty(refindx),
    ft_error('the requested reference channel is not found in the data');
  end
  if any(strncmp(tok, 'rpt', 3)),
    Cr = freq.crsspctrm(:,chanindx,refindx);
    Pr = freq.crsspctrm(:,refindx,refindx);
  else
    Cr = freq.crsspctrm(chanindx,refindx);
    Pr = freq.crsspctrm(refindx,refindx);
  end
end  

% do a sanity check on the cross-spectral-density matrix
if any(isnan(Cf(:)))
  ft_error('The cross-spectral-density matrix is not complete');
end
if any(isnan(Cr(:)))
  ft_error('The cross-spectral-density with the reference channel is not complete');
end

