function [freq] = ft_freqinterpolate(cfg, freq)

% FT_FREQINTERPOLATE interpolates frequencies by looking at neighbouring
% values or simply replaces a piece in the spectrum by NaN.
%
% Use as
%   freq = ft_freqinterpolate(cfg, freq)
% where freq is the output of FT_FREQANALYSIS or FT_FREQDESCRIPTIVES and the
% configuration may contain
%   cfg.method   = 'nan', 'linear' (default = 'nan')
%   cfg.foilim   = Nx2 matrix with begin and end of each interval to be
%                  interpolated (default = [49 51; 99 101; 149 151])

% Undocumented local options:
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%   cfg.outputfile = one can specify output as file to save to disk


% Copyright (C) 2009, Aldemar Torres Valderama
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

fieldtripdefs

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% set the default values
if ~isfield(cfg, 'method'),     cfg.method = 'nan';                     end
if ~isfield(cfg, 'foilim'),     cfg.foilim = [49 51; 99 101; 149 151];  end
if ~isfield(cfg, 'inputfile'),  cfg.inputfile                   = [];    end
if ~isfield(cfg, 'outputfile'), cfg.outputfile                  = [];    end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    freq = loadvar(cfg.inputfile, 'freq');
  end
end

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'datatype', 'freq', 'feedback', 'yes');

for i = 1:size(cfg.foilim,1)
  % determine the exact frequency bins to interpolate
  peakbeg = nearest(freq.freq, cfg.foilim(i,1));
  peakend = nearest(freq.freq, cfg.foilim(i,2));
  % update the configuration
  cfg.foilim(i,1) = freq.freq(peakbeg);
  cfg.foilim(i,2) = freq.freq(peakend);
  
  if strcmp(cfg.method, 'nan')
    switch freq.dimord
      case 'chan_freq'
        freq.powspctrm(:,peakbeg:peakend) = nan;
      case 'chan_freq_time'
        freq.powspctrm(:,peakbeg:peakend,:) = nan;
      case 'rpt_chan_freq'
        freq.powspctrm(:,:,peakbeg:peakend) = nan;
      case 'rpt_chan_freq_time'
        freq.powspctrm(:,:,peakbeg:peakend,:) = nan;
      otherwise
        error('unsupported dimord');
    end % switch
    
  elseif strcmp(cfg.method, 'linear')
    error('not yet implemented');
    
  else
    error('unsupported method');
  end
end % for each frequency range

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = freq.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'freq', freq); % use the variable name "data" in the output file
end
