function type = ft_chantype(input, desired)

% FT_CHANTYPE determines for each individual channel what type of data it
% represents, e.g. a planar gradiometer, axial gradiometer, magnetometer,
% trigger channel, etc. If you want to know what the acquisition system is
% (e.g. ctf151 or neuromag306), you should not use this function but
% FT_SENSTYPE instead.
%
% Use as
%   type = ft_chantype(hdr)
%   type = ft_chantype(sens)
%   type = ft_chantype(label)
% or as
%   type = ft_chantype(hdr,   desired)
%   type = ft_chantype(sens,  desired)
%   type = ft_chantype(label, desired)
%
% If the desired unit is not specified as second input argument, this
% function returns a Nchan*1 cell-array with a string describing the type
% of each channel.
%
% If the desired unit is specified as second input argument, this function
% returns a Nchan*1 boolean vector with "true" for the channels of the
% desired type and "false" for the ones that do not match.
%
% The specification of the channel types depends on the acquisition system,
% for example the ctf275 system includes the following type of channels:
% meggrad, refmag, refgrad, adc, trigger, eeg, headloc, headloc_gof.
%
% See also FT_READ_HEADER, FT_SENSTYPE, FT_CHANUNIT

% Copyright (C) 2008-2015, Robert Oostenveld
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

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

% this is to avoid a recursion loop
persistent recursion
if isempty(recursion)
  recursion = false;
end

if nargin<2
  desired = [];
end

% determine the type of input, this is handled similarly as in FT_CHANUNIT
isheader = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'Fs');
isgrad   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  &&  isfield(input, 'ori'); % old style
iselec   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  && ~isfield(input, 'ori'); % old style
isgrad   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'coilpos')) || isgrad;             % new style 
iselec   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'elecpos')) || iselec;             % new style 
islabel  = isa(input, 'cell')    && ~isempty(input) && isa(input{1}, 'char');

if isheader
  % this speeds up the caching in real-time applications
  input.nSamples = 0;
end

current_argin = {input, desired};
if isequal(current_argin, previous_argin)
  % don't do the type detection again, but return the previous output from cache
  type = previous_argout{1};
  return
end

if isheader
  label = input.label;
  numchan = length(label);
elseif isgrad
  label   = input.label;
  numchan = length(label);
elseif iselec
  label   = input.label;
  numchan = length(label);
elseif islabel
  label   = input;
  numchan = length(label);
elseif isfield(input, 'label')
  % this is a last resort: I don't know what it is, but perhaps the labels are informative
  label   = input.label;
  numchan = length(label);
else
  error('the input that was provided to this function cannot be deciphered');
end

if isfield(input, 'chantype')
  % start with the provided channel types
  type = input.chantype(:);
else
  % start with unknown type for all channels
  type = repmat({'unknown'}, numchan, 1);
end

if ft_senstype(input, 'unknown')
  % don't bother doing all subsequent checks to determine the type of sensor array
  
elseif isheader && (ft_senstype(input, 'neuromag') || ft_senstype(input, 'babysquid74')) 
  % channames-KI is the channel kind, 1=meg, 202=eog, 2=eeg, 3=trigger (I am not sure, but have inferred this from a single test file)
  % chaninfo-TY is the Coil type (0=magnetometer, 1=planar gradiometer)
  if isfield(input, 'orig') && isfield(input.orig, 'channames')
    for sel=find(input.orig.channames.KI(:)==202)'
      type{sel} = 'eog';
    end
    for sel=find(input.orig.channames.KI(:)==2)'
      type{sel} = 'eeg';
    end
    for sel=find(input.orig.channames.KI(:)==3)'
      type{sel} = 'digital trigger';
    end
    % determine the MEG channel subtype
    selmeg=find(input.orig.channames.KI(:)==1)';
    for i=1:length(selmeg)
      if input.orig.chaninfo.TY(i)==0
        type{selmeg(i)} = 'megmag';
      elseif input.orig.chaninfo.TY(i)==1
        % FIXME this might also be a axial gradiometer in case the BabySQUID data is read with the old reading routines
        type{selmeg(i)} = 'megplanar';
      end
    end
    
  elseif isfield(input, 'orig') && isfield(input.orig, 'chs') && isfield(input.orig.chs, 'coil_type')
    % all the chs.kinds and chs.coil_types are obtained from the MNE manual, p.210-211
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==2)' % planar gradiometers
      type(sel) = {'megplanar'}; %Neuromag-122 planar gradiometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3012)' %planar gradiometers
      type(sel) = {'megplanar'}; %Type T1 planar grad
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3013)' %planar gradiometers
      type(sel) = {'megplanar'}; %Type T2 planar grad
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3014)' %planar gradiometers
      type(sel) = {'megplanar'}; %Type T3 planar grad
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3022)' %magnetometers
      type(sel) = {'megmag'};    %Type T1 magenetometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3023)' %magnetometers
      type(sel) = {'megmag'};    %Type T2 magenetometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3024)' %magnetometers
      type(sel) = {'megmag'};    %Type T3 magenetometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==7001)' %axial gradiometer
      type(sel) = {'megaxial'};
    end
    for sel=find([input.orig.chs.kind]==301)' %MEG reference channel, located far from head
      type(sel) = {'ref'};
    end
    for sel=find([input.orig.chs.kind]==2)'   %EEG channels
      type(sel) = {'eeg'};
    end
    for sel=find([input.orig.chs.kind]==201)' %MCG channels
      type(sel) = {'mcg'};
    end
    for sel=find([input.orig.chs.kind]==3)' %Stim channels
      if any([input.orig.chs(sel).logno] == 101) %new systems: 101 (and 102, if enabled) are digital; low numbers are 'pseudo-analog' (if enabled)
        type(sel([input.orig.chs(sel).logno] == 101)) = {'digital trigger'};
        type(sel([input.orig.chs(sel).logno] == 102)) = {'digital trigger'};
        type(sel([input.orig.chs(sel).logno] <= 32))  = {'analog trigger'};
        others = [input.orig.chs(sel).logno] > 32 & [input.orig.chs(sel).logno] ~= 101 & ...
          [input.orig.chs(sel).logno] ~= 102;
        type(sel(others)) = {'other trigger'};
      elseif any([input.orig.chs(sel).logno] == 14) %older systems: STI 014/015/016 are digital; lower numbers 'pseudo-analog'(if enabled)
        type(sel([input.orig.chs(sel).logno] == 14)) = {'digital trigger'};
        type(sel([input.orig.chs(sel).logno] == 15)) = {'digital trigger'};
        type(sel([input.orig.chs(sel).logno] == 16)) = {'digital trigger'};
        type(sel([input.orig.chs(sel).logno] <= 13)) = {'analog trigger'};
        others = [input.orig.chs(sel).logno] > 16;
        type(sel(others)) = {'other trigger'};
      else
        warning('There does not seem to be a suitable trigger channel.');
        type(sel) = {'other trigger'};
      end
    end
    for sel=find([input.orig.chs.kind]==202)' %EOG
      type(sel) = {'eog'};
    end
    for sel=find([input.orig.chs.kind]==302)' %EMG
      type(sel) = {'emg'};
    end
    for sel=find([input.orig.chs.kind]==402)' %ECG
      type(sel) = {'ecg'};
    end
    for sel=find([input.orig.chs.kind]==502)' %MISC
      type(sel) = {'misc'};
    end
    for sel=find([input.orig.chs.kind]==602)' %Resp
      type(sel) = {'respiration'};
    end
  end
  
elseif ft_senstype(input, 'babysquid74')
  % the name can be something like "MEG 001" or "MEG001" or "MEG 0113" or "MEG0113"
  % i.e. with two or three digits and with or without a space
  sel = myregexp('^MEG', label);
  type(sel) = {'megaxial'};
  
elseif ft_senstype(input, 'neuromag122')
  % the name can be something like "MEG 001" or "MEG001" or "MEG 0113" or "MEG0113"
  % i.e. with two or three digits and with or without a space
  sel = myregexp('^MEG', label);
  type(sel) = {'megplanar'};
  
elseif ft_senstype(input, 'neuromag306') && isgrad
  % there should be 204 planar gradiometers and 102 axial magnetometers
  if isfield(input, 'tra')
    tmp = sum(abs(input.tra),2);
    sel = (tmp==median(tmp));
    type(sel) = {'megplanar'};
    sel = (tmp~=median(tmp));
    type(sel) = {'megmag'};
  end
  
elseif ft_senstype(input, 'ctf') && isheader
  % According to one source of information meg channels are 5, refmag 0,
  % refgrad 1, adcs 18, trigger 11, eeg 9.
  %
  % According to another source of information it is as follows
  %     refMagnetometers: 0
  %      refGradiometers: 1
  %             meg_sens: 5
  %             eeg_sens: 9
  %                  adc: 10
  %             stim_ref: 11
  %           video_time: 12
  %                  sam: 15
  %     virtual_channels: 16
  %             sclk_ref: 17
  
  % start with an empty one
  origSensType = [];
  if isfield(input, 'orig')
    if isfield(input.orig, 'sensType') && isfield(input.orig, 'Chan')
      % the header was read using the open-source MATLAB code that originates from CTF and that was modified by the FCDC
      origSensType = input.orig.sensType;
    elseif isfield(input.orig, 'res4') && isfield(input.orig.res4, 'senres')
      % the header was read using the CTF p-files, i.e. readCTFds
      origSensType =  [input.orig.res4.senres.sensorTypeIndex];
    elseif isfield(input.orig, 'sensor') && isfield(input.orig.sensor, 'info')
      % the header was read using the CTF importer from the NIH and Daren Weber
      origSensType = [input.orig.sensor.info.index];
    end
  end
  
  if isempty(origSensType)
    warning('could not determine channel type from the CTF header');
  end
  
  for sel=find(origSensType(:)==5)'
    type{sel} = 'meggrad';
  end
  for sel=find(origSensType(:)==0)'
    type{sel} = 'refmag';
  end
  for sel=find(origSensType(:)==1)'
    type{sel} = 'refgrad';
  end
  for sel=find(origSensType(:)==18)'
    type{sel} = 'adc';
  end
  for sel=find(origSensType(:)==11)'
    type{sel} = 'trigger';
  end
  for sel=find(origSensType(:)==17)'
    type{sel} = 'clock';
  end
  for sel=find(origSensType(:)==9)'
    type{sel} = 'eeg';
  end
  for sel=find(origSensType(:)==29)'
    type{sel} = 'reserved'; % these are "reserved for future use", but relate to head localization
  end
  for sel=find(origSensType(:)==13)'
    type{sel} = 'headloc'; % these represent the x, y, z position of the head coils
  end
  for sel=find(origSensType(:)==28)'
    type{sel} = 'headloc_gof'; % these represent the goodness of fit for the head coils
  end
  % for sel=find(origSensType(:)==23)'
  %   type{sel} = 'SPLxxxx'; % I have no idea what these are
  % end
  
elseif ft_senstype(input, 'ctf') && isgrad
  % in principle it is possible to look at the number of coils, but here the channels are identified based on their name
  sel = myregexp('^M[ZLR][A-Z][0-9][0-9]$', input.label);
  type(sel) = {'meggrad'};            % normal gradiometer channels
  sel = myregexp('^S[LR][0-9][0-9]$', input.label);
  type(sel) = {'meggrad'};            % normal gradiometer channels in the 64 channel CTF system
  sel = myregexp('^B[GPQR][0-9]$', input.label);
  type(sel) = {'refmag'};             % reference magnetometers
  sel = myregexp('^[GPQR][0-9][0-9]$', input.label);
  type(sel) = {'refgrad'};            % reference gradiometers
  
elseif ft_senstype(input, 'ctf') && islabel
  % the channels have to be identified based on their name alone
  sel = myregexp('^M[ZLR][A-Z][0-9][0-9]$', label);
  type(sel) = {'meggrad'};                % normal gradiometer channels
  sel = myregexp('^S[LR][0-9][0-9]$', label);
  type(sel) = {'meggrad'};            % normal gradiometer channels in the 64 channel CTF system
  sel = myregexp('^B[GPR][0-9]$', label);
  type(sel) = {'refmag'};             % reference magnetometers
  sel = myregexp('^[GPQR][0-9][0-9]$', label);
  type(sel) = {'refgrad'};            % reference gradiometers
  
elseif ft_senstype(input, 'bti')
  if isfield(input, 'orig') && isfield(input.orig, 'config')
    configname = {input.orig.config.channel_data.name};
    configtype = [input.orig.config.channel_data.type];
    
    if ~isequal(configname(:), input.label(:))
      % reorder the channels according to the order in input.label
      [sel1, sel2] = match_str(input.label, configname);
      configname = configname(sel2);
      configtype = configtype(sel2);
      configdata = input.orig.config.channel_data(sel2);
    end
    numloops = zeros(size(configdata));
    for i=1:length(configdata)
      if isfield(configdata(i).device_data, 'total_loops')
        numloops(i) = configdata(i).device_data.total_loops;
      end
    end
    
    % these are taken from bti2grad
    type(configtype==1 & numloops==1) = {'megmag'};
    type(configtype==1 & numloops==2) = {'meggrad'};
    type(configtype==2) = {'eeg'};
    type(configtype==3) = {'ref'}; % not known if mag or grad
    type(configtype==4) = {'aux'};
    type(configtype==5) = {'trigger'};
    
    % refine the distinction between refmag and refgrad to make the types
    % in grad and header consistent
    sel = myregexp('^M[CLR][xyz][aA]*$', label);
    type(sel) = {'refmag'};
    sel = myregexp('^G[xyz][xyz]A$', label);
    type(sel) = {'refgrad'};    
  else
    % determine the type on the basis of the channel labels
    % all 4D-BTi MEG channels start with "A" followed by a number
    % all 4D-BTi reference channels start with M or G
    % all 4D-BTi EEG channels start with E, except for the 248-MEG/32-EEG system in Warsaw where they end with -1
    sel = myregexp('^A[0-9]+$', label);
    type(sel) = {'meg'};
    sel = myregexp('^M[CLR][xyz][aA]*$', label);
    type(sel) = {'refmag'};
    sel = myregexp('^G[xyz][xyz]A$', label);
    type(sel) = {'refgrad'};
    
    if isgrad && isfield(input, 'tra')
      gradtype = repmat({'unknown'}, size(input.label));
      gradtype(strncmp('A', input.label, 1)) = {'meg'};
      gradtype(strncmp('M', input.label, 1)) = {'refmag'};
      gradtype(strncmp('G', input.label, 1)) = {'refgrad'};
      % look at the number of coils of the meg channels
      selchan = find(strcmp('meg', gradtype));
      for k = 1:length(selchan)
        ncoils = length(find(input.tra(selchan(k),:)==1));
        if ncoils==1,
          gradtype{selchan(k)} = 'megmag';
        elseif ncoils==2,
          gradtype{selchan(k)} = 'meggrad';
        end
      end
      [selchan, selgrad] = match_str(label, input.label);
      type(selchan) = gradtype(selgrad);
    end
    
    % deal with additional channel types based on the names
    if isheader && issubfield(input, 'orig.channel_data.chan_label')
      tmplabel = {input.orig.channel_data.chan_label};
      tmplabel = tmplabel(:);
    else
      tmplabel = label; % might work
    end
    sel = find(strcmp('unknown', type));
    if ~isempty(sel)
      type(sel) = ft_chantype(tmplabel(sel));
      sel       = find(strcmp('unknown', type));
      if ~isempty(sel)
        % channels that start with E are assumed to be EEG
        % channels that end with -1 are also assumed to be EEG, see http://bugzilla.fcdonders.nl/show_bug.cgi?id=2389
        type(sel(cellfun(@(x) strcmp(x(end-1:end),'-1') || strcmp(x(1),'E'), label(sel)))) = {'eeg'};
      end
    end
  end
  
elseif ft_senstype(input, 'itab') && isheader
  origtype = [input.orig.ch.type];
  type(origtype==0) = {'unknown'};
  type(origtype==1) = {'ele'};
  type(origtype==2) = {'mag'}; % might be magnetometer or gradiometer, look at the number of coils
  type(origtype==4) = {'ele ref'};
  type(origtype==8) = {'mag ref'};
  type(origtype==16) = {'aux'};
  type(origtype==32) = {'param'};
  type(origtype==64) = {'digit'};
  type(origtype==128) = {'flag'};
  % these are the channels that are visible to fieldtrip
  chansel = 1:input.orig.nchan;
  type = type(chansel);
  
elseif ft_senstype(input, 'yokogawa') && isheader
  % This is to recognize Yokogawa channel types from the original header
  % This is from the original documentation
  NullChannel                   = 0;
  MagnetoMeter                  = 1;
  AxialGradioMeter              = 2;
  PlannerGradioMeter            = 3;
  RefferenceChannelMark         = hex2dec('0100');
  RefferenceMagnetoMeter        = bitor( RefferenceChannelMark, MagnetoMeter      );
  RefferenceAxialGradioMeter    = bitor( RefferenceChannelMark, AxialGradioMeter  );
  RefferencePlannerGradioMeter  = bitor( RefferenceChannelMark, PlannerGradioMeter);
  TriggerChannel                = -1;
  EegChannel                    = -2;
  EcgChannel                    = -3;
  EtcChannel                    = -4;
  if ft_hastoolbox('yokogawa_meg_reader')
    % shorten names
    ch_info = input.orig.channel_info.channel;
    type_orig = [ch_info.type];
    
    sel = (type_orig == NullChannel);
    type(sel) = {'null'};
    sel = (type_orig == MagnetoMeter);
    type(sel) = {'megmag'};
    sel = (type_orig == AxialGradioMeter);
    type(sel) = {'meggrad'};
    sel = (type_orig == PlannerGradioMeter);
    type(sel) = {'megplanar'};
    sel = (type_orig == RefferenceMagnetoMeter);
    type(sel) = {'refmag'};
    sel = (type_orig == RefferenceAxialGradioMeter);
    type(sel) = {'refgrad'};
    sel = (type_orig == RefferencePlannerGradioMeter);
    type(sel) = {'refplanar'};
    sel = (type_orig == TriggerChannel);
    type(sel) = {'trigger'};
    sel = (type_orig == EegChannel);
    type(sel) = {'eeg'};
    sel = (type_orig == EcgChannel);
    type(sel) = {'ecg'};
    sel = (type_orig == EtcChannel);
    type(sel) = {'etc'};
  elseif ft_hastoolbox('yokogawa')
    sel = (input.orig.channel_info(:, 2) == NullChannel);
    type(sel) = {'null'};
    sel = (input.orig.channel_info(:, 2) == MagnetoMeter);
    type(sel) = {'megmag'};
    sel = (input.orig.channel_info(:, 2) == AxialGradioMeter);
    type(sel) = {'meggrad'};
    sel = (input.orig.channel_info(:, 2) == PlannerGradioMeter);
    type(sel) = {'megplanar'};
    sel = (input.orig.channel_info(:, 2) == RefferenceMagnetoMeter);
    type(sel) = {'refmag'};
    sel = (input.orig.channel_info(:, 2) == RefferenceAxialGradioMeter);
    type(sel) = {'refgrad'};
    sel = (input.orig.channel_info(:, 2) == RefferencePlannerGradioMeter);
    type(sel) = {'refplanar'};
    sel = (input.orig.channel_info(:, 2) == TriggerChannel);
    type(sel) = {'trigger'};
    sel = (input.orig.channel_info(:, 2) == EegChannel);
    type(sel) = {'eeg'};
    sel = (input.orig.channel_info(:, 2) == EcgChannel);
    type(sel) = {'ecg'};
    sel = (input.orig.channel_info(:, 2) == EtcChannel);
    type(sel) = {'etc'};
  end
  
elseif ft_senstype(input, 'yokogawa') && isgrad
  % all channels in the gradiometer definition are meg
  % type(1:end) = {'meg'};
  % channels are identified based on their name: only magnetic as isgrad==1
  sel = myregexp('^M[0-9][0-9][0-9]$', input.label);
  type(sel) = {'megmag'};
  sel = myregexp('^AG[0-9][0-9][0-9]$', input.label);
  type(sel) = {'meggrad'};
  sel = myregexp('^PG[0-9][0-9][0-9]$', input.label);
  type(sel) = {'megplanar'};
  sel = myregexp('^RM[0-9][0-9][0-9]$', input.label);
  type(sel) = {'refmag'};
  sel = myregexp('^RAG[0-9][0-9][0-9]$', input.label);
  type(sel) = {'refgrad'};
  sel = myregexp('^RPG[0-9][0-9][0-9]$', input.label);
  type(sel) = {'refplanar'};
  
elseif ft_senstype(input, 'yokogawa') && islabel
  % the yokogawa channel labels are a mess, so autodetection is not possible
  % type(1:end) = {'meg'};
  sel = myregexp('[0-9][0-9][0-9]$', label);
  type(sel) = {'null'};
  sel = myregexp('^M[0-9][0-9][0-9]$', label);
  type(sel) = {'megmag'};
  sel = myregexp('^AG[0-9][0-9][0-9]$', label);
  type(sel) = {'meggrad'};
  sel = myregexp('^PG[0-9][0-9][0-9]$', label);
  type(sel) = {'megplanar'};
  sel = myregexp('^RM[0-9][0-9][0-9]$', label);
  type(sel) = {'refmag'};
  sel = myregexp('^RAG[0-9][0-9][0-9]$', label);
  type(sel) = {'refgrad'};
  sel = myregexp('^RPG[0-9][0-9][0-9]$', label);
  type(sel) = {'refplanar'};
  sel = myregexp('^TRIG[0-9][0-9][0-9]$', label);
  type(sel) = {'trigger'};
  sel = myregexp('^EEG[0-9][0-9][0-9]$', label);
  type(sel) = {'eeg'};
  sel = myregexp('^ECG[0-9][0-9][0-9]$', label);
  type(sel) = {'ecg'};
  sel = myregexp('^ETC[0-9][0-9][0-9]$', label);
  type(sel) = {'etc'};
  
elseif ft_senstype(input, 'itab') && isheader
  sel = ([input.orig.ch.type]==0);
  type(sel) = {'unknown'};
  sel = ([input.orig.ch.type]==1);
  type(sel) = {'unknown'};
  sel = ([input.orig.ch.type]==2);
  type(sel) = {'megmag'};
  sel = ([input.orig.ch.type]==8);
  type(sel) = {'megref'};
  sel = ([input.orig.ch.type]==16);
  type(sel) = {'aux'};
  sel = ([input.orig.ch.type]==64);
  type(sel) = {'digital'};
  % not all channels are actually processed by fieldtrip, so only return
  % the types fopr the ones that read_header and read_data return
  type = type(input.orig.chansel);
  
elseif ft_senstype(input, 'itab') && isgrad
  % the channels have to be identified based on their name alone
  sel = myregexp('^MAG_[0-9][0-9][0-9]$', label);
  type(sel) = {'megmag'};
  sel = myregexp('^MAG_[0-9][0-9]$', label); % for the itab28 system
  type(sel) = {'megmag'};
  sel = myregexp('^MAG_[0-9]$', label); % for the itab28 system
  type(sel) = {'megmag'};
  sel = myregexp('^REF_[0-9][0-9][0-9]$', label);
  type(sel) = {'megref'};
  sel = myregexp('^AUX.*$', label);
  type(sel) = {'aux'};
  
elseif ft_senstype(input, 'itab') && islabel
  % the channels have to be identified based on their name alone
  sel = myregexp('^MAG_[0-9][0-9][0-9]$', label);
  type(sel) = {'megmag'};
  sel = myregexp('^REF_[0-9][0-9][0-9]$', label);
  type(sel) = {'megref'};
  sel = myregexp('^AUX.*$', label);
  type(sel) = {'aux'};
  
elseif ft_senstype(input, 'eeg') && islabel
  % use an external helper function to define the list with EEG channel names
  type(match_str(label, ft_senslabel('eeg1005'))) = {'eeg'};          % this includes all channels from the 1010 and 1020 arrangement
  type(match_str(label, ft_senslabel(ft_senstype(input)))) = {'eeg'}; % this will work for biosemi, egi and other detected channel arrangements

elseif ft_senstype(input, 'eeg') && iselec
  % all channels in an electrode definition must be eeg channels
  type(:) = {'eeg'};
    
elseif ft_senstype(input, 'eeg') && isheader
  % use an external helper function to define the list with EEG channel names
  type(match_str(input.label, ft_senslabel(ft_senstype(input)))) = {'eeg'};
  
elseif ft_senstype(input, 'plexon') && isheader
  % this is a complete header that was read from a Plexon *.nex file using read_plexon_nex
  for i=1:numchan
    switch input.orig.VarHeader(i).Type
      case 0
        type{i} = 'spike';
      case 1
        type{i} = 'event';
      case 2
        type{i} = 'interval';  % Interval variables?
      case 3
        type{i} = 'waveform';
      case 4
        type{i} = 'population'; % Population variables ?
      case 5
        type{i} = 'analog';
      otherwise
        % keep the default 'unknown' type
    end
  end
  
end % ft_senstype

% if possible, set additional types based on channel labels
label2type = {
  {'ecg', 'ekg'};
  {'emg'};
  {'eog', 'heog', 'veog'};
  {'lfp'};
  {'eeg'};
  {'trigger', 'trig', 'dtrig'};
  };
for i = 1:numel(label2type)
  for j = 1:numel(label2type{i})
    type(intersect(strmatch(label2type{i}{j}, lower(label)), find(strcmp(type, 'unknown')))) = label2type{i}(1);
  end
end

if all(strcmp(type, 'unknown')) && ~recursion
  % try whether only lowercase channel labels makes a difference
  if islabel
    recursion = true;
    type = ft_chantype(lower(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = lower(input.label);
    recursion = true;
    type = ft_chantype(input);
    recursion = false;
  end
end

if all(strcmp(type, 'unknown')) && ~recursion
  % try whether only uppercase channel labels makes a difference
  if islabel
    recursion = true;
    type = ft_chantype(upper(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = upper(input.label);
    recursion = true;
    type = ft_chantype(input);
    recursion = false;
  end
end

if nargin>1
  % return a boolean vector
  if isequal(desired, 'meg') || isequal(desired, 'ref')
    % only compare the first three characters, i.e. meggrad or megmag should match
    type = strncmp(desired, type, 3);
  else
    type = strcmp(desired, type);
  end
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {type};
previous_argin  = current_argin;
previous_argout = current_argout;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function match = myregexp(pat, list)
match = false(size(list));
for i=1:numel(list)
  match(i) = ~isempty(regexp(list{i}, pat, 'once'));
end
