function [dataout] = ft_electrodermalactivity(cfg, datain)

% FT_ELECTRODERMALACTIVITY estimates the electrodermal activity from a recording of
% the electric resistance of the skin.
%
% Use as
%   eda = ft_electrodermalactivity(cfg, data)
% where the input data is a structure as obtained from FT_PREPROCESSING.
%
% The configuration structure has the following options
%   cfg.channel        = selected channel for processing, see FT_CHANNELSELECTION
%   cfg.feedback       = 'yes' or 'no'
%   cfg.medianwindow   = scalar, length of window for median filter in seconds (default = 8)
%
% After using this function you can use FT_REDEFINETRIAL and FT_TIMELOCKANLAYSIS to
% investigate electrodermal responses (EDRs) to stimulation. You can use
% FT_ARTIFACT_THRESHOLD to determine the timing and frequency of nonspecific EDRs.
%
% See https://doi.org/10.1111/j.1469-8986.2012.01384.x "Publication recommendations
% for electrodermal measurements" by the SPR for an introduction in electrodermal
% methods and for recommendations.
%
% See also FT_HEARTRATE, FT_HEADMOVEMENT, FT_REGRESSCONFOUND

% Copyright (C) 2018, Robert Oostenveld, DCCN
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    datain
ft_preamble provenance datain

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% check if the input data is valid for this function, the input data must be raw
datain = ft_checkdata(datain, 'datatype', 'raw', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729

% set the default options
cfg.channel        = ft_getopt(cfg, 'channel', {});
cfg.feedback       = ft_getopt(cfg, 'feedback', 'yes');
cfg.medianwindow   = ft_getopt(cfg, 'medianwindow', 8);   % in seconds
cfg.preproc        = ft_getopt(cfg, 'preproc', []);

% copy some of the fields over to the new data structure
dataout = keepfields(datain, {'time', 'fsample', 'sampleinfo', 'trialinfo'});
dataout.label = {'tonic', 'phasic'};
dataout.trial = {};  % this is to be determined in the main code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.channel = ft_channelselection(cfg.channel, datain.label);
assert(numel(cfg.channel)==1, 'you should specify exactly one channel');

fsample = datain.fsample;
chansel = strcmp(datain.label, cfg.channel{1});
medianwindow = round(cfg.medianwindow*fsample); % in samples

for trllop=1:numel(datain.trial)
  dat   = datain.trial{trllop}(chansel,:);
  label = datain.label(chansel);
  time  = datain.time{trllop};

  if ~isempty(cfg.preproc)
    % apply the preprocessing to the selected channel
    [dat, label, time, cfg.preproc] = preproc(dat, label, time, cfg.preproc, 0, 0);
  end

  tonic  = ft_preproc_medianfilter(dat, medianwindow);
  phasic = dat - tonic;

  if istrue(cfg.feedback)
    figure
    subplot(3,1,1)
    plot(time, dat)
    title('preprocessed')
    subplot(3,1,2)
    plot(time, tonic)
    title('tonic')
    subplot(3,1,3)
    plot(time, phasic)
    title('phasic')
  end

  dataout.trial{trllop} = [tonic; phasic];
end % for trllop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout
