function fieldtrip2fiff(filename, data)

% FIELDTRIP2FIFF saves a FieldTrip raw data structure as a fiff-file, allowing it
% to be further analyzed by the Neuromag/Elekta/Megin software, or in MNE-python.
%
% Use as
%   fieldtrip2fiff(filename, data)
% where filename is the name of the output file, and data is a raw data structure
% as obtained from FT_PREPROCESSING, or a timelock structure obtained from
% FT_TIMELOCKANALYSIS.
%
% If the data comes from preprocessing and has only one trial, then it writes the
% data into raw continuous format. If present in the data, the original header
% is reused (also removing the non-used channels). Otherwise, the function 
% attempts to create the header, which might or might not be correct (e.g. with
% respect to the scaling and the sensor locations. If the data contains events in 
% the cfg structure, it writes the events in the MNE format (three columns) into 
% a file based on "filename", ending with "-eve.fif"
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
fifffile  = fullfile(pathstr ,[name '.fif']);
eventfile = fullfile(pathstr ,[name '-eve.fif']);

% ensure the mne-toolbox to be on the path
ft_hastoolbox('mne', 1);

% check if the input data is valid for this function
data   = ft_checkdata(data, 'datatype', {'raw', 'timelock'}, 'hassampleinfo', 'yes', 'feedback', 'yes');
istlck = ft_datatype(data, 'timelock');
isepch = ft_datatype(data, 'raw');
israw  = false;
if isepch && numel(data.trial) == 1
  isepch = false;
  israw  = true;
end
if istlck
  fsample = 1./mean(diff(data.time));
else
  fsample = 1./mean(diff(data.time{1}));
end

% Create a fiff-header, or take information from the original header if possible
if isfield(data, 'hdr') && isfield(data.hdr, 'orig') && isfield(data.hdr.orig, 'meas_id')
  ft_notice('The data contains header information that seems to be derived from a FIFF file,\nre-using header information, channel locations may be read \nfrom .grad and .elec in the data');
  info = data.hdr.orig;
else
  info.meas_id.version = NaN;
  info.meas_id.machid  = [NaN;NaN];
  info.meas_id.secs    = NaN;
  info.meas_id.usecs   = NaN;
  info.meas_date       = [NaN;NaN];
  info.acq_pars        = []; % needed by raw
  info.acq_stim        = []; % needed by raw
  info.highpass        = NaN;
  info.lowpass         = NaN;
  
  % these are not strictly necessary, but the inverse functions in MNE works better if this matrix is present
  info.dev_head_t.from  = 1;
  info.dev_head_t.to    = 4;
  info.dev_head_t.trans = eye(4);
  
  info.ctf_head_t = [];
  info.dig        = [];
  info.projs      = struct('kind', {}, 'data', {}, 'active', {}, 'desc', {});
  info.comps      = struct('kind', {}, 'data', {}, 'ctfkind', {}, 'save_calibrated', {}, 'rowcals', {}, 'colcals', {});
  info.bads       = [];
  
  info.isaverage    = double(istlck);
  info.isepoched    = double(isepch);
  info.iscontinuous = double(israw);
  info.sfreq        = fsample;
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
  
  ft_warning('writing epochs to MNE is not implemented yet');
  
  ntrl    = numel(data.trial);
  aspect_kind = num2cell(ones(ntrl,1)*100);
  is_smsh = num2cell(zeros(ntrl,1)); % FIXME: How could we tell? + Don't know what this is
  nave    = num2cell(ones(ntrl,1));  % FIXME: Use the real value
  first   = num2cell(data.sampleinfo(:, 1));
  last    = num2cell(data.sampleinfo(:, 2));
  comment = cell(ntrl,1);
  for j = 1:ntrl
    comment{j,1} = sprintf('FieldTrip data, category/trial %d', j);
  end
  evoked = struct('aspect_kind', aspect_kind, 'is_smsh', is_smsh, 'nave', nave, ...
    'first', first, 'last', last, 'comment', comment, 'times', data.time(:), 'epochs', data.trial(:));

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
if isfield(data, 'hdr') && isfield(data.hdr, 'orig') && isfield(data.hdr.orig, 'chs')
  
  ft_warning('Using the original channel information from the header, this may be wrong if channels have been pruned, reordered, rescaled');
  [i_label, i_chs] = match_str(data.label, {data.hdr.orig.chs.ch_name}');
  chs(i_label)     = data.hdr.orig.chs(i_chs);
else

  % otherwise reconstruct a usable channel information array
  ft_warning('Reconstructing channel information, based on the available data, it might be inaccurate\n');
  FIFF = fiff_define_constants; % some constants are not defined in the MATLAB function
  
  hasgrad = isfield(data, 'grad');
  haselec = isfield(data, 'elec');
  nchan   = numel(data.label);
  
  stype = zeros(nchan, 1).*-1;
  indx  = nan(nchan,1);
  
  if haselec
    [i_labeeg, i_elec] = match_str(data.label, data.elec.label);
    stype(i_labeeg) = 0; %meg/eeg/other = 1/0/-1
    indx(i_labeeg)  = i_elec; %indexes into the grad/elec
   
    data.elec = ft_convert_units(data.elec, 'cm'); % MNE uses cm for EEG (FIXME is this true)
  end
  
  if hasgrad
    [i_labmeg, i_grad] = match_str(data.label, data.grad.label);
    stype(i_labmeg) = 1;
    indx(i_labmeg)  = i_grad;
  
    data.grad = undobalancing(data.grad); % if this fails, then it's difficult to write out meaningful channel info to begin with
    coiltype  = grad2coiltype(data.grad);
    coilunit  = grad2coilunit(data.grad, FIFF);
  end

  scanno  = num2cell((1:nchan)');
  ch_name = data.label(:);
  range   = num2cell(ones(nchan,1));
  cal     = num2cell(ones(nchan,1));

  
  cnt_grad = 0;
  cnt_elec = 0;
  cnt_else = 0;
  chs = struct('scanno', scanno, 'ch_name', ch_name, 'range', range, 'cal', cal)';
  for k = 1:nchan
  
    switch stype(k)
      case 1
        % MEG
        cnt_grad = cnt_grad + 1;
        
        chs(1,k).kind      = FIFF.FIFFV_MEG_CH;
        chs(1,k).logno     = cnt_grad;
        chs(1,k).coil_type = coiltype(indx(k));
        chs(1,k).unit      = coilunit(indx(k));
        
        pos = data.grad.chanpos(indx(k),:);
        ori = data.grad.chanori(indx(k),:);

        this_z = ori;          
        this_z = this_z / norm(this_z);
        this_x = cross([0 0 1], this_z);
        if all(this_x==0)
          this_x = [1 0 0];
        else
          this_x = this_x / norm(this_x);
        end
        this_y = cross(this_z, this_x);
        R      = [this_x(:) this_y(:) this_z(:)];

        chs(1,k).coil_trans   = [R pos(:); 0 0 0 1];
        chs(1,k).unit_mul     = 0;
        chs(1,k).coord_frame  = FIFF.FIFFV_COORD_HEAD;
        chs(1,k).eeg_loc      = [];
        chs(1,k).loc          = [pos(:); R(:)];
      case 0
        % EEG
        cnt_elec = cnt_elec + 1;
       
        chs(1,k).kind         = FIFF.FIFFV_EEG_CH;
        chs(1,k).logno        = cnt_elec;
        chs(1,k).coil_type    = NaN;
        chs(1,k).coil_trans   = [];
        chs(1,k).unit         = 107; % volts FIFF.FIFF_UNIT_V
        chs(1,k).unit_mul     = -6; % micro FIFF.FIFF_UNITM_MU
        chs(1,k).coord_frame  = FIFF.FIFFV_COORD_DEVICE;
        chs(1,k).eeg_loc      = [elec.chanpos(i_elec,:)' zeros(3,1)] / 100;
        chs(1,k).loc          = [chs(1,k).eeg_loc(:); 0; 1; 0; 0; 0; 1];

      case -1
        % OTHER
        cnt_else = cnt_else + 1;
        
        chs(1,k).kind         = NaN;
        chs(1,k).logno        = cnt_else;
        chs(1,k).coil_type    = NaN;
        chs(1,k).coil_trans   = [];
        chs(1,k).unit         = NaN;
        chs(1,k).unit_mul     = 0;
        chs(1,k).coord_frame  = NaN;
        chs(1,k).eeg_loc      = [];
        chs(1,k).loc          = zeros(12,1);

      otherwise
    end

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

function coiltype = grad2coiltype(grad)

stype = ft_senstype(grad);
def   = mne_load_coil_def('coil_def.dat');
def   = def([def.accuracy]==0);

coiltype = [];
switch stype
  case 'neuromag122'
  case 'neuromag306'
  case {'ctf151' 'ctf275'}
  case 'bti148'
  case 'bti248'
  otherwise
    
end

% treat as point magnetometer system
sel = strcmp({def.description}', 'Point magnetometer');
if isempty(coiltype)
  coiltype = ones(numel(grad.label),1).*def(sel).id;
end

function coilunit = grad2coilunit(grad, FIFF)

coilunit = zeros(numel(grad.label),1)-1;
for k = 1:numel(grad.chanunit)
  switch grad.chanunit{k}
    case 'T'
      coilunit(k) = FIFF.FIFF_UNIT_T;
    case 'T/m'
      coilunit(k) = FIFF.FIFF_UNIT_T_M;
    otherwise
      % figure out what to do if T/cm or so -> I think that this goes into
      % the unit_mul field
  end
end

