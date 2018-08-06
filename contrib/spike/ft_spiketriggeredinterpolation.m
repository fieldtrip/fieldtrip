function [data] = ft_spiketriggeredinterpolation(cfg, data)

% FT_SPIKETRIGGEREDINTERPOLATION interpolates the data in the LFP channels
% around the spikes that are detected in the spike channels, or replaces
% the LFP around the spike with NaNs. The purpose of this procedure is to
% allow analysis of spikes and LFPs recorded from the same electrode, as
% the spike energy would bleed in the LFP.
%
% Use as
%   [data] = ft_spiketriggeredinterpolation(cfg, data)
%
% The input data should be organised in a structure as obtained from the
% FT_PREPROCESSING function. The configuration should be according to
%
%   cfg.method       = string, The interpolation method can be 'nan',
%                     'cubic', 'linear', 'nearest', spline', 'pchip'
%                     (default = 'nan'). See INTERP1 for more details.
%   cfg.timwin       = [begin end], duration of LFP segment around each spike (default =
%                      [-0.005 0.005]) to be removed
%   cfg.interptoi    = value, time in seconds used for interpolation, which
%                      must be larger than timwin (default = 0.2)
%   cfg.spikechannel = string, name of single spike channel to trigger on
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')
%
% The output will contain all channels of the input, only the data in the
% selected channels will be interpolated or replaced with NaNs.
%
% See also FT_SPIKETRIGGEREDSPECTRUM, FT_SPIKETRIGGEREDAVERAGE

% Copyright (C) 2008, Thilo Womelsdorf
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
ft_preamble provenance data
ft_preamble trackconfig

% check input data structure
data = ft_checkdata(data,'datatype', 'raw', 'feedback', 'yes');

% these were supported in the past, but are not any more (for consistency with other spike functions)
cfg = ft_checkconfig(cfg, 'forbidden', {'inputfile', 'outputfile'});  

%get the options
cfg.timwin         = ft_getopt(cfg, 'timwin',[-0.005 0.005]);
cfg.spikechannel   = ft_getopt(cfg,'spikechannel', []);
cfg.channel        = ft_getopt(cfg,'channel', 'all');
cfg.feedback       = ft_getopt(cfg,'feedback', 'yes');
cfg.method         = ft_getopt(cfg,'method', 'nan');
cfg.interptoi      = ft_getopt(cfg,'interptoi',0.2);

% ensure that the options are valid
cfg = ft_checkopt(cfg,'timwin','ascendingdoublebivector');
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'channel', {'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'feedback', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'method', 'char');
cfg = ft_checkopt(cfg,'interptoi', 'doublescalar');

if strcmp(cfg.method, 'nan'), cfg.interptoi = 0; end

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=1:nchans
    spikechan(j) = spikechan(j) + all(data.trial{i}(j,:)==0 | data.trial{i}(j,:)==1 | data.trial{i}(j,:)==2);
  end
end
spikechan = (spikechan==ntrial);

% determine the channels for interpolation
cfg.channel = ft_channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);

% determine the spike channel on which will be triggered
cfg.spikechannel = ft_channelselection(cfg.spikechannel, data.label);
spikesel         = match_str(data.label, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel);    % number of spike channels

if nspikesel==0, error('no spike channel selected'); end
if nspikesel>1, error('only supported for a single spike channel'); end

if ~spikechan(spikesel)
  error('the selected spike channel seems to contain continuous data');
end

% this determines the segment that will be replaced around each spike
begpad = round(cfg.timwin(1)*data.fsample);
endpad = round(cfg.timwin(2)*data.fsample);

% this determines the segment that is used for the inperpolation around each spike
interppad = round( cfg.interptoi*data.fsample);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate the LFP around the spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:ntrial
  spikesmp = find(data.trial{i}(spikesel,:));
  
  fprintf('Spike Triggered Interpolation: processing trial %d of %d (%d spikes)\n', i, ntrial, length(spikesmp));
  
  ft_progress('init', cfg.feedback, 'interpolating spikes');
  for j=1:length(spikesmp)
    ft_progress(i/ntrial, 'interpolating spike %d of %d\n', j, length(spikesmp));
    begsmp = spikesmp(j) + begpad;
    endsmp = spikesmp(j) + endpad;
    
    begsmp_interp = begsmp - interppad;
    endsmp_interp = endsmp + interppad;
    
    if begsmp_interp<1
      continue,
    end
    if endsmp_interp>size(data.trial{i},2)
      continue,
    end
    
    if strcmp(cfg.method,'nan')
      % only replace with NaNs
      data.trial{i}(chansel,begsmp:endsmp) = NaN;
      
    else
      % interpolate the data around the spike
      xall  = [begsmp_interp          : endsmp_interp];
      x     = [begsmp_interp:begsmp-1   endsmp+1:endsmp_interp];
      y     =  data.trial{i}(chansel,x) ;
      if length(chansel) > 1; y = y'; end %our channels are in rows but interp1 need cols
      yi    = interp1(x,y,xall,cfg.method);
      if length(chansel) > 1; yi = yi'; end %get us back to rows
      
      % store the interpolated segment back in the data
      data.trial{i}(chansel,xall) = yi;            
    end % if strcmp(cfg.method)
    
  end % for each spike in this trial
  ft_progress('close');
  
end % for each trial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the data has been modified on the fly, only update the configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data

