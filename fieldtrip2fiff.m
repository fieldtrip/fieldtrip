function fieldtrip2fiff(filename, data)

% FIELDTRIP2FIFF saves a fieldtrip data structure as a fiff-file, in order
% to be useable by the Neuromag software, or in MNE suite
%
% Use as
%  fieldtrip2fiff(filename, data)
% 
% where filename is the name of the output file, and data is a fieldtrip
% raw data structure, or a timelock structure.
%

% Copyright (C) 2012 Jan-Mathijs Schoffelen
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
% $Id: ft_fieldtrip2fiff.m 5187 2012-01-31 08:42:56Z jansch $

revision = '$Id: ft_fieldtrip2fiff.m 5187 2012-01-31 08:42:56Z jansch $';

ft_defaults                 % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble help            % this will show the function help if nargin==0 and return an error
ft_preamble callinfo        % this records the time and memory usage at teh beginning of the function

% ensure that the filename has the correct extension
[pathstr,name,ext] = fileparts(filename);
if isempty(ext),
  filename = [filename, '.fif'];
elseif ~strcmp(ext, 'fif')
  error('if the filename is specified with extension, this should read .fif');
end

% ensure the mne-toolbox to be on the path
ft_hastoolbox('mne', 1);

% check the input data
data   = ft_checkdata(data, 'datatype', {'raw', 'timelock'}, 'feedback', 'yes');
istlck = ft_datatype(data, 'timelock');
israw  = ft_datatype(data, 'raw');

% Create a fiff-header, or take it from the original header if possible
if ft_senstype(data, 'neuromag') && isfield(data, 'hdr')
  % Tune the original FIFF header to match with the data
  info       = data.hdr.orig;
  info.sfreq = fsample;
  if info.nchan ~= size(dataout{1},1);
    fprintf('WARNING: A non-matching number of channels in the FT data structure and original file header. Integrity of the data might be compromised\');
    info.nchan = size(dataout{1},1);
    info.chs = info.chs(1:size(dataout{1},1)); % FIXME: Terrible hack to tolerate removal of channels
  end
else
  info.meas_id.version = nan;
  info.meas_id.machid  = [nan;nan];
  info.meas_id.secs    = nan;
  info.meas_id.usecs   = nan;
  info.meas_date       = [nan;nan];
  
  info.nchan    = numel(data.label);
  info.highpass = nan;
  info.lowpass  = nan;
  info.dev_head_t = [];
  info.ctf_head_t = [];
  info.dig      = [];
  info.projs    = [];
  info.comps    = [];
  info.bads     = [];
  info.ch_names = data.label(:)';
  info.chs      = grad2fiff(data);
  if istlck,
    info.sfreq     = 1./mean(diff(data.time));
    info.isaverage = 1;
    info.isepoched = 0;
    info.iscontinuous = 0;
  elseif israw,
    info.sfreq     = 1./mean(diff(data.time{1}));
    info.isaverage = 0;
    info.isepoched = 1;
    info.iscontinuous = 0;
  end
end

if istlck
  evoked.aspect_kind = 100;
  evoked.is_smsh     = 0;
  evoked.nave        = max(data.dof(:));
  evoked.first       = round(data.time(1)*info.sfreq);
  evoked.last        = round(data.time(end)*info.sfreq);
  evoked.times       = data.time;
  evoked.comment     = sprintf('FieldTrip data averaged');
  evoked.epochs      = data.avg;
elseif israw
  for j = 1:length(dataout)
    evoked(j).aspect_kind = 100;
    evoked(j).is_smsh     = 0; % FIXME: How could we tell?
    evoked(j).nave        = 1; % FIXME: Use the real value
    evoked(j).first       = round(data.time{j}(1)*info.sfreq);
    evoked(j).last        = round(data.time{j}(end)*info.sfreq);
    evoked(j).times       = data.time{j};
    evoked(j).comment     = sprintf('FieldTrip data, category/trial %d', j);
    evoked(j).epochs      = data.trial{j};
  end  
end

fiffdata.info   = info;
fiffdata.evoked = evoked;

fiff_write_evoked(filename, fiffdata);

ft_postamble callinfo         % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and matlab version etc. to the output cfg
ft_postamble previous datain  % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble history dataout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar dataout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option

%-------------------
% subfunction
function [chs] = grad2fiff(data)

if isfield(data, 'grad')
  hasgrad = true;
else
  hasgrad = false;
end

for k = 1:numel(data.label)
  % create a struct for each channel
  chs(1,k).scanno       = 1;
  chs(1,k).logno        = 1;
  chs(1,k).kind         = nan;
  chs(1,k).range        = 1;
  chs(1,k).cal          = 1;
  chs(1,k).coil_type    = nan;
  chs(1,k).loc          = zeros(12,1);
  chs(1,k).coil_trans   = eye(4);
  chs(1,k).eeg_loc      = [];
  chs(1,k).coord_frame  = nan;
  chs(1,k).unit         = nan;
  chs(1,k).unit_mul     = 0;
  chs(1,k).ch_name      = data.label{k};
end

if hasgrad
  % the position information can be recovered, too
  [i1,i2] = match_str(data.label, data.grad.label);
  for k = 1:numel(i1)
    chs(1,i1(k)).kind        = 1;
    chs(1,i1(k)).coil_type   = 3022;
    chs(1,i1(k)).unit        = 112;
    chs(1,i1(k)).coord_frame = int32(4); % head coordinates
    chs(1,i1(k)).loc(1:3)    = data.grad.chanpos(i2(k),:);
    chs(1,i1(k)).loc(4:end)  = reshape(eye(3),[9 1]);
  end
end

% Below the documentation to the original code contributed by Lauri
% Parkkonen
%
% fieldtrip2fiff(filename, data)
%
% Write the data in the FieldTrip structure 'data' to a FIFF file
% 'filename'. If an average
%
% Caveats:
%
% - the FIFF tag indicating the number of trials in average is set to unity
%     as there is no generic way to determine a proper value.
%
% - the FIFF tag indicating the use of 'MaxShield' is set to 'no'.
%
% - if channels have been removed in FieldTrip, channel information for
%     a matching number of channels is copied from the original FIFF
%     header. If channels have been removed from anywhere else except from
%     the end, this simple hack screws up the pairing of signals and channel
%     identities. BEWARE!

% (C)opyright Lauri Parkkonen, 2010 - 2011
%
% $Log$
%

% % Construct the evoked data aspect for each category
% ave = isfield(data, 'average');
% trl = isfield(data, 'trial');
%
% if ave && trl
%     fprintf('WARNING: Data structure contains both an average and individual trials; writing out the average');
%     dataout = data.average;
% elseif ave && ~trl
%     dataout = data.average;
% elseif ~ave && trl
%     dataout = data.trial;
% else
%     error('No average or trial data to write out');
% end
% 
% % Tune the original FIFF header to match with the data
% info = data.hdr.orig;
% info.sfreq = data.fsample;
% if info.nchan ~= size(dataout{1},1);
%     fprintf('WARNING: A non-matching number of channels in the FT data structure and original file header. Integrity of the data might be compromised\');
%     info.nchan = size(dataout{1},1);
%     info.chs = info.chs(1:size(dataout{1},1)); % FIXME: Terrible hack to tolerate removal of channels
% end
% 
% 
% for c = 1:length(dataout)
%     evoked(c).aspect_kind = 100;
%     evoked(c).is_smsh = 0; % FIXME: How could we tell?
%     evoked(c).nave = 1; % FIXME: Use the real value
%     evoked(c).first = round(data.time{1}(1) * data.fsample);
%     evoked(c).last = round(data.time{1}(end) * data.fsample);
%     evoked(c).times = data.time{1};
%     evoked(c).comment = sprintf('FieldTrip data, category/trial %d', c);
%     evoked(c).epochs = dataout{c};
% end
% 
% fiffdata.info = info;
% fiffdata.evoked = evoked;
% 
% fiff_write_evoked(filename, fiffdata);
% 
% end
