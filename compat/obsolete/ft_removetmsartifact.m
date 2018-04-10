function data = ft_removetmsartifact(cfg, data)

% FT_REMOVETMSARTIFACT removes TMS artifacts from EEG data
%
% %%
% NOTE: Please be aware that this function is deprecated. Please follow the
% TMS-EEG tutorial instead at http://www.fieldtriptoolbox.org/tutorial/tms-eeg
% %%
%
% Use as
%  data = ft_removetmsartifact(cfg, data)
% where the input data is a raw data, for example obtained from FT_PREPROCESSING, and
% cfg is a configuratioun structure that should contain
%   cfg.method      = string, can be 'twopassfilter', 'interpolatepulse'
%   cfg.pulseonset  = value or vector, time in seconds of the TMS pulse in seconds
%
% The following options pertain to the 'replace' method
%   cfg.pulsewidth  = value, pulse pulsewidth to be removed in seconds
%   cfg.offset      = value, offset with respect to pulse onset to start
%                     replacing, in seconds.
%
% The following options pertain to the 'twopassfilter' method
%   cfg.lpfreq      = number in Hz
%   cfg.lpfiltord   = lowpass  filter order
%   cfg.lpfilttype  = digital filter type, 'but' or 'fir' or 'firls' (default = 'but')
%   cfg.pulsewidth  = value, pulse pulsewidth to be removed in seconds. If
%                     set to 0, entire trial will be filtered.
%   cfg.offset      = value, offset with respect to pulse onset to start
%                     filtering, in seconds.
%
% See also FT_REJECTARTIFACT, FT_REJECTCOMPONENT

% Copyrights (C) 2012, Robert Oostenveld
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

% DEPRECATED by jimher on 19 September 2013
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1791 for more details
warning('FT_REMOVETMSARTIFACT is deprecated, please follow TMS-EEG tutorial instead (http://www.fieldtriptoolbox.org/tutorial/tms-eeg).')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar datain
ft_preamble provenance datain
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw'}, 'feedback', 'yes');

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'method'});

% get the options
cfg.method     = ft_getopt(cfg, 'method');        % there is no default
cfg.pulseonset = ft_getopt(cfg, 'pulseonset');
cfg.pulsewidth = ft_getopt(cfg, 'pulsewidth');
cfg.lpfiltord  = ft_getopt(cfg, 'lpfiltord', 2);
cfg.lpfilttype = ft_getopt(cfg, 'lpfilttype', 'but');
cfg.lpfreq     = ft_getopt(cfg, 'lpfreq', 30);
cfg.offset     = ft_getopt(cfg, 'offset', 0);
cfg.fillmethod = ft_getopt(cfg, 'fillmethod');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numtrl  = length(data.trial);
temp_pulse = [];

if ~isfield(data, 'fsample')
  fsample = 1/mean(diff(data.time{1}));
else
  fsample = data.fsample;
end

if isnumeric(cfg.pulsewidth) && numel(cfg.pulsewidth)==1; temp_pulse = cfg.pulsewidth; end

% copy for all trials
if isnumeric(cfg.pulseonset) && numel(cfg.pulseonset)==1; cfg.pulseonset = repmat(cfg.pulseonset, numtrl, 1); end
if isnumeric(cfg.pulsewidth) && numel(cfg.pulsewidth)==1; cfg.pulsewidth = repmat(cfg.pulsewidth, numtrl, 1); end

% check wether fields are cell where necessary
if ~iscell(cfg.pulseonset); cfg.pulseonset = num2cell(cfg.pulseonset); end
if ~iscell(cfg.pulsewidth); cfg.pulsewidth = num2cell(cfg.pulsewidth); end

if isempty(cfg.pulseonset) || isempty(cfg.pulsewidth)
  for i=1:numtrl
    [onset, width] = pulsedetect(data.trial{i});
    % these should be expressed in seconds
    cfg.pulseonset{i} = data.time{i}(onset);

    if ~isempty(temp_pulse)
      cfg.pulsewidth{i} = repmat(temp_pulse, 1, length(onset));
    else
      cfg.pulsewidth{i} = width;
    end

    fprintf('detected %d pulses in trial %d\n', length(onset), i);
  end
end % estimate pulse onset and width

switch cfg.method

  case 'twopassfilter'
    for i=1:numtrl
      for j=1:length(cfg.pulseonset{i})
        %tmssample = nearest(data.time{i}, cfg.pulseonset{i}(j));
        pulseonset = cfg.pulseonset{i}(j);
        pulsewidth = cfg.pulsewidth{i}(j);
        offset = cfg.offset;

        % express it in samples,
        pulseonset = nearest(data.time{i}, pulseonset);
        pulsewidth  = round(pulsewidth*fsample);
        offset = round(offset*fsample);

        % get the part of the data that is left and right of the TMS pulse artifact
        dat1 = data.trial{i}(:,1:pulseonset);
        dat2 = data.trial{i}(:,(pulseonset+1:end));

        % filter the two pieces, prevent filter artifacts
        [filt1] = ft_preproc_lowpassfilter(dat1,fsample,cfg.lpfreq,cfg.lpfiltord,cfg.lpfilttype,'onepass');
        [filt2] = ft_preproc_lowpassfilter(dat2,fsample,cfg.lpfreq,cfg.lpfiltord,cfg.lpfilttype,'onepass-reverse');

        % stitch the left and right parts of the data back together
        %data.trial{i} = [filt1 filt2];
        fill = [filt1 filt2];

        % determine a short window around the tms pulse
        begsample = pulseonset + offset;
        endsample = pulseonset + pulsewidth + offset - 1;

        % replace data in the pulse window with a filtered version
        if pulsewidth == 0
          data.trial{i} = fill;
        else
          data.trial{i}(:,begsample:endsample) = fill(:,begsample:endsample);
        end
      end % for pulses
    end % for trials

  case 'interpolatepulse'
    for i=1:numtrl
      for j=1:length(cfg.pulseonset{i})

        pulseonset = cfg.pulseonset{i}(j);
        pulsewidth = cfg.pulsewidth{i}(j);
        offset = cfg.offset;

        % express it in samples,
        pulseonset = nearest(data.time{i}, pulseonset);
        pulsewidth  = round(pulsewidth*fsample);
        offset = round(offset*fsample);

        begsample = pulseonset + offset;
        endsample = pulseonset + pulsewidth + offset - 1;

        % determine a short window before the TMS pulse
        begsample1 = begsample - pulsewidth;
        endsample1 = begsample - 1;

        % determine a short window after the TMS pulse
        begsample2 = endsample + 1;
        endsample2 = endsample + pulsewidth;

        dat1 = data.trial{i}(:,begsample1:endsample1);
        dat2 = data.trial{i}(:,begsample2:endsample2);
        %fill = dat1(:,randperm(size(dat1,2))); % randomly shuffle the data points
        %fill = mean(dat1,2) + cumsum(std(dat1,[],2).*randn(size(dat1,1),size(dat1,2)));
%         fill = linspace(mean(dat1,2),mean(dat2,2),endsample1-begsample1+1);
%         fill = fill + cumsum(std(dat1,[],2).*randn(size(dat1,1),size(dat1,2)));

%         fill = cumsum(std(dat1,[],2).*randn(size(dat1,1),size(dat1,2)));


        switch cfg.fillmethod
            case 'fft'
            fft_dat1 = fft(dat1);
            fft_dat2 = fft(dat2);
            fill = real(ifft(mean([fft_dat1; fft_dat2])));
    %         fill = std(dat1,[],2).*randn(size(dat1,1),size(dat1,2));
    %         fill = fill .* repmat(hann(size(fill,2))',size(fill,1),1);
    %         %fill = fill - linspace(fill(:,1),fill(:,end),endsample1-begsample1+1);
    % %         fill = fill + linspace(mean(dat1,2),mean(dat2,2),endsample1-begsample1+1);
    %         fill = fill + linspace(dat1(:,end),dat2(:,1),endsample1-begsample1+1);
            case 'zeros'
                fill = zeros(size(dat1,1),size(dat1,2));
            case 'randperm'
                fill = dat1(:,randperm(size(dat1,2))); % randomly shuffle the data points
            case 'brown'
                fill = linspace(mean(dat1,2),mean(dat2,2),endsample1-begsample1+1);
                fill = fill + cumsum(std(dat1,[],2).*randn(size(dat1,1),size(dat1,2)));
            case 'linear'
                fill = interp1([1:size(dat1,2) 2*size(dat1,2)+1:3*size(dat1,2)], [dat1 dat2]', size(dat1,2)+1:2*size(dat1,2),'linear')';
            case 'linear+noise'
                fill = interp1([1:size(dat1,2) 2*size(dat1,2)+1:3*size(dat1,2)], [dat1 dat2]', size(dat1,2)+1:2*size(dat1,2),'linear')';
                fill = fill(2:end-1) + std(dat1,[],2).*randn(size(dat1,1),size(dat1,2));
            case 'spline'
                fill = interp1([1:size(dat1,2) 2*size(dat1,2)+1:3*size(dat1,2)], [dat1 dat2]', size(dat1,2)+1:2*size(dat1,2),'spline')';
            case 'cubic'
                fill = interp1([1:size(dat1,2) 2*size(dat1,2)+1:3*size(dat1,2)], [dat1 dat2]', size(dat1,2)+1:2*size(dat1,2),'cubic')';
            case 'cubic+noise'
                fill = interp1([1:size(dat1,2) 2*size(dat1,2)+1:3*size(dat1,2)], [dat1 dat2]', size(dat1,2)+1:2*length(dat1),'cubic')';
                fill = fill + std(dat1,[],2).*randn(size(dat1,1),size(dat1,2));
        end

            % FIXME an alternative would be to replace it with an interpolated version of the signal just around it
            % FIXME an alternative would be to replace it with nan
            % FIXME an alternative would be to replace it with random noise

            % replace the data in the pulse window with a random shuffled version of the data just around it
        data.trial{i}(:,begsample:endsample) = fill;

      end % for pulses
    end % for trials

  otherwise
    ft_error('unsupported method');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that detects the onset and pulsewidth of one or multiple TMS pulses
% that are present as artifact in a segment of multi-channel EEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [onset, pulsewidth] = pulsedetect(dat)
[nchan, ntime] = size(dat);
for i=1:nchan
  dat(i,:) = dat(i,:) - median(dat(i,:));
end
dat = sum(abs(dat),1);
threshold = 0.5 * max(dat);
dat = dat>threshold;
dat = [0 diff(dat) 0];
onset  = find(dat== 1);
offset = find(dat==-1) - 1;
pulsewidth  = offset - onset;
% make the pulse a bit wider
offset = offset - 2*pulsewidth;
pulsewidth = pulsewidth*5;
