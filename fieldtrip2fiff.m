function fieldtrip2fiff(filename, data)

% FIELDTRIP2FIFF saves a FieldTrip raw data structure as a fiff-file, allowing it
% to be further analyzed by the Elekta/Neuromag software, or in the MNE suite
% software.
%
% Use as
%   fieldtrip2fiff(filename, data)
% where filename is the name of the output file, and data is a raw data structure
% as obtained from FT_PREPROCESSING, or a timelock structure obtained from
% FT_TIMELOCKANALYSIS.
%
% If the data comes from preprocessing and has only one trial, then it writes the
% data into raw continuous format. If present in the data, the original header
% from neuromag is reused (also removing the non-used channels). Otherwise, the
% function tries to create a correct header, which might or might not contain the
% correct scaling and channel location. If the data contains events in the cfg
% structure, it writes the events in the MNE format (three columns) into a file
% based on "filename", ending with "-eve.fif"
%
% See also FT_DATATYPE_RAW, FT_DATATYPE_TIMELOCK

% Copyright (C) 2012-2013, Jan-Mathijs Schoffelen, Gio Piantoni
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

% this ensures that the path is correct and that the ft_defaults global variable is available
ft_defaults

% ensure that the filename has the correct extension
[pathstr, name, ext] = fileparts(filename);
if ~strcmp(ext, '.fif')
  ft_error('if the filename is specified with extension, this should read .fif');
end
fifffile = fullfile(pathstr ,[name '.fif']);
eventfile = fullfile(pathstr ,[name '-eve.fif']);

% ensure the mne-toolbox to be on the path
ft_hastoolbox('mne', 1);

% check if the input data is valid for this function
data   = ft_checkdata(data, 'datatype', {'raw', 'timelock'}, 'feedback', 'yes');
istlck = ft_datatype(data, 'timelock');
isepch = ft_datatype(data, 'raw');
israw = false;
if isepch && numel(data.trial) == 1
  isepch = false;
  israw = true;
end

% Create a fiff-header, or take it from the original header if possible
if ft_senstype(data, 'neuromag') && isfield(data, 'hdr')
  fprintf('Using the original FIFF header, but channel locations are read \nfrom .grad and .elec in the data, if they exist\n')
  info = data.hdr.orig;
  
else
  info.meas_id.version = NaN;
  info.meas_id.machid  = [NaN;NaN];
  info.meas_id.secs    = NaN;
  info.meas_id.usecs   = NaN;
  info.meas_date       = [NaN;NaN];
  info.acq_pars = []; % needed by raw
  info.acq_stim = []; % needed by raw
  
  info.highpass = NaN;
  info.lowpass  = NaN;
  % no strictly necessary, but the inverse functions in MNE works better if
  % this matrix is present
  info.dev_head_t.from = 1;
  info.dev_head_t.to = 4;
  info.dev_head_t.trans = eye(4);
  
  info.ctf_head_t = [];
  info.dig      = [];
  info.projs = struct('kind', {}, 'active', {}, 'desc', {}, 'data', {});
  info.comps = struct('ctfkind', {}, 'kind', {}, 'save_calibrated', {}, ...
    'rowcals', {}, 'colcals', {}, 'data', {});
  info.bads     = [];
  
  if isepch
    info.sfreq     = 1./mean(diff(data.time{1}));
    info.isaverage = 0;
    info.isepoched = 1;
    info.iscontinuous = 0;
    
  elseif istlck
    info.sfreq     = 1./mean(diff(data.time));
    info.isaverage = 1;
    info.isepoched = 0;
    info.iscontinuous = 0;
    
  end
  
end

if israw
  info.sfreq = data.fsample;
elseif isepch
  info.sfreq = 1 ./ mean(diff(data.time{1}));
elseif istlck
  info.sfreq = 1 ./ mean(diff(data.time));
end

info.ch_names = data.label(:)';
info.chs      = sens2fiff(data);
info.nchan    = numel(data.label);

if israw
  [outfid, cals] = fiff_start_writing_raw(fifffile, info);
  fiff_write_raw_buffer(outfid, data.trial{1}, cals);
  fiff_finish_writing_raw(outfid);
  
  % write events, if they exists
  if isfield(data, 'cfg')
    event = ft_findcfg(data.cfg, 'event');
  else
    event = [];
  end
  if ~isempty(event)
    eve = convertevent(event);
    mne_write_events(eventfile, eve);
    fprintf('Writing events to %s\n', eventfile)
  end
  
elseif isepch
  
  ft_error('writing epochs to MNE is not implemented yet')
  
  for j = 1:length(data.trial)
    evoked(j).aspect_kind = 100;
    evoked(j).is_smsh     = 0; % FIXME: How could we tell?
    evoked(j).nave        = 1; % FIXME: Use the real value
    evoked(j).first       = round(data.time{j}(1)*info.sfreq);
    evoked(j).last        = round(data.time{j}(end)*info.sfreq);
    evoked(j).times       = data.time{j};
    evoked(j).comment     = sprintf('FieldTrip data, category/trial %d', j);
    evoked(j).epochs      = data.trial{j};
  end
  
  % fiffdata.info   = info;
  % fiffdata.evoked = evoked;
  % fiff_write_XXX(fifffile, fiffdata);
  
elseif istlck
  evoked.aspect_kind = 100;
  evoked.is_smsh     = 0;
  evoked.nave        = max(data.dof(:));
  evoked.first       = round(data.time(1)*info.sfreq);
  evoked.last        = round(data.time(end)*info.sfreq);
  evoked.times       = data.time;
  evoked.comment     = sprintf('FieldTrip data averaged');
  evoked.epochs      = data.avg;
  
  fiffdata.info   = info;
  fiffdata.evoked = evoked;
  fiff_write_evoked(fifffile, fiffdata);
  
end

%-------------------
% subfunction
function [chs] = sens2fiff(data)

% use orig information if available
if isfield(data, 'hdr') && isfield(data.hdr, 'orig') && ...
    isfield(data.hdr.orig, 'chs')
  
  [dummy, i_label, i_chs] = intersect(data.label, {data.hdr.orig.chs.ch_name});
  chs(i_label) = data.hdr.orig.chs(i_chs);
  return
  
end

% otherwise reconstruct it
fprintf('Reconstructing channel locations, it might be inaccurate\n')
FIFF = fiff_define_constants; % some constants are not defined in the MATLAB function

if isfield(data, 'grad')
  hasgrad = true;
else
  hasgrad = false;
end

if isfield(data, 'elec')
  haselec = true;
  elec = ft_convert_units(data.elec, 'cm'); % that MNE uses cm
else
  haselec = false;
end

cnt_grad = 0;
cnt_elec = 0;
cnt_else = 0;

for k = 1:numel(data.label)
  % create a struct for each channel
  chs(1,k).scanno       = k;
  chs(1,k).ch_name      = data.label{k};
  
  chs(1,k).range        = 1;
  chs(1,k).cal          = 1;
  
  i_grad = false;
  i_elec = false;
  if hasgrad
    i_grad = strcmp(data.grad.label, data.label{k});
  elseif haselec
    i_elec = strcmp(elec.label, data.label{k});
  end
  
  if any(i_grad)
    chs(1,k).kind         = FIFF.FIFFV_MEG_CH;
    cnt_grad = cnt_grad + 1;
    chs(1,k).logno        = cnt_grad;
    switch data.grad.chantype{i_grad}
      case 'megmag'
        chs(1,k).coil_type    = 3024;
        chs(1,k).unit         = FIFF.FIFF_UNIT_T;
      case 'megplanar'
        chs(1,k).coil_type    = 3012;
        chs(1,k).unit         = FIFF.FIFF_UNIT_T_M;
      case 'meggrad'
        chs(1,k).coil_type    = 3022;
        chs(1,k).unit         = FIFF.FIFF_UNIT_T;
      otherwise
        fprintf('Unknown channel type %s, assigned to meggrad', data.grad.chantype{i_grad})
        chs(1,k).coil_type    = 3022;
        chs(1,k).unit         = FIFF.FIFF_UNIT_T;
    end
    chs(1,k).coil_trans   = eye(4);
    chs(1,k).unit_mul     = 0;
    chs(1,k).coord_frame  = FIFF.FIFFV_COORD_HEAD;
    chs(1,k).eeg_loc      = [];
    chs(1,k).loc          = [data.grad.chanpos(i_grad,:)'; reshape(eye(3),[9 1])];
    
  elseif any(i_elec)
    chs(1,k).kind         = FIFF.FIFFV_EEG_CH;
    cnt_elec = cnt_elec + 1;
    chs(1,k).logno        = cnt_elec;
    chs(1,k).coil_type    = NaN;
    chs(1,k).coil_trans   = [];
    chs(1,k).unit         = 107; % volts FIFF.FIFF_UNIT_V
    chs(1,k).unit_mul     = -6; % micro FIFF.FIFF_UNITM_MU
    chs(1,k).coord_frame  = FIFF.FIFFV_COORD_DEVICE;
    chs(1,k).eeg_loc      = [elec.chanpos(i_elec,:)' zeros(3,1)] / 100;
    chs(1,k).loc          = [chs(1,k).eeg_loc(:); 0; 1; 0; 0; 0; 1];
    
  else
    chs(1,k).kind         = NaN;
    cnt_else = cnt_else + 1;
    chs(1,k).logno        = cnt_else;
    chs(1,k).coil_type    = NaN;
    chs(1,k).coil_trans   = [];
    chs(1,k).unit         = NaN;
    chs(1,k).unit_mul     = 0;
    chs(1,k).coord_frame  = NaN;
    chs(1,k).eeg_loc      = [];
    chs(1,k).loc          = zeros(12,1);
    
  end
  
end

function eve = convertevent(event)
% tentative code, with lots of assumption

%CTF should use backpanel trigger
backpanel = strcmp({event.type}, 'backpanel trigger');
if any(backpanel)
  fprintf('Writing the value of the backpanel trigger into the event file\n')
  trigger = [event(backpanel).value];
  
  eve = zeros(numel(trigger), 3);
  eve(:,1) = [event(backpanel).sample];
  eve(:,3) = [event(backpanel).value];
  return
end

% use ev_type and ev_value
ev_type = unique({event.type});
% convert to cell of strings
if any(cellfun(@isnumeric, {event.value}))
  event_value = cellfun(@num2str, {event.value}, 'uni', false);
else
  event_value = {event.value};
end
ev_value = unique(event_value);

eve = zeros(numel(event), 3);

for i1 = 1:numel(ev_type)
  for i2 = 1:numel(ev_value)
    i_type = strcmp({event.type}, ev_type{i1});
    i_value = strcmp(event_value, ev_value{i2});
    % if events are numeric & there's only one event type keep original code:
    if ~isempty(str2num(ev_value{i2})) && numel(ev_type) == 1
        marker = str2num(ev_value{i2});
    else
      marker = i1 * 10 + i2;
    end
    
    if any(i_type & i_value)
      eve(i_type & i_value, 1) = [event(i_type & i_value).sample];
      eve(i_type & i_value, 3) = marker;
    end
    
  end
end

% report event coding
newev = unique(eve(:,3));
if all(cellfun(@isnumeric, {event.value})) && numel(ev_type) == 1
    fprintf('EVENT codes remain the same.\n')
else
  fprintf('EVENTS have been coded as:\n')
  for i = 1:numel(newev)
    i_type = floor(newev(i)/10);
    i_value = mod(newev(i), 10);
    fprintf('type: %s, value %s -> % 3d\n', ev_type{i_type}, ev_value{i_value}, newev(i))
  end
end
