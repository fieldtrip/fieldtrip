function [freq] = freqinterpolate(cfg, freq)

% FREQINTERPOLATE interpolates frequencies by looking at neighbouring
% values or simply replaces a piece in the spectrum by NaN.
%
% Use as
%   freq = freqinterpolate(cfg, freq)
% where freq is the output of FREQANALYSIS or FREQDESCRIPTIVES and the
% configuration may contain
%   cfg.method   = 'nan', 'linear' (default = 'nan')
%   cfg.foilim   = Nx2 matrix with begin and end of each interval to be
%                  interpolated (default = [49 51; 99 101; 149 151])

% Copyright (C) 2009, Aldemar Torres Valderama
%
% $Log: freqinterpolate.m,v $
% Revision 1.1  2009/02/24 09:24:32  roboos
% first implementation of new functino, can be used for plotting of powerspectra that contain line-noise and harmonic peaks
%

fieldtripdefs

% check if the input data is valid for this function
freq = checkdata(freq, 'datatype', 'freq', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the default values
if ~isfield(cfg, 'method'),     cfg.method = 'nan';                     end
if ~isfield(cfg, 'foilim'),     cfg.foilim = [49 51; 99 101; 149 151];  end

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

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: freqinterpolate.m,v 1.1 2009/02/24 09:24:32 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = freq.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;
