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
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_FREQANALYSIS, FT_FREQDESCRIPTIVES, FT_FREQSIMULATION

% Copyright (C) 2009, Aldemar Torres Valderama
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
ft_preamble loadvar freq
ft_preamble provenance freq
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'datatype', 'freq', 'feedback', 'yes');

% set the defaults
cfg.method = ft_getopt(cfg, 'method', 'nan');
cfg.foilim = ft_getopt(cfg, 'foilim', [49 51; 99 101; 149 151]);


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
        ft_error('unsupported dimord');
    end % switch

  elseif strcmp(cfg.method, 'linear')
    ft_error('not yet implemented');

  else
    ft_error('unsupported method');
  end
end % for each frequency range

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   freq
ft_postamble provenance freq
ft_postamble history    freq
ft_postamble savevar    freq
