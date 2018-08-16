function chantype = ft_chantype(input, desired)

% FT_CHANTYPE determines for each individual channel what chantype of data it
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
% for example the ctf275 system includes the following tyoe of channels:
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
isdata   = isa(input, 'struct')  && ~isheader && (isfield(input, 'hdr') || isfield(input, 'grad') || isfield(input, 'elec') || isfield(input, 'opto'));
isgrad   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  &&  isfield(input, 'ori'); % old style
iselec   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  && ~isfield(input, 'ori'); % old style
isgrad   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'coilpos')) || isgrad;             % new style
iselec   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'elecpos')) || iselec;             % new style
isopto   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'transceiver');
islabel  = isa(input, 'cell')    && ~isempty(input) && isa(input{1}, 'char');

if isheader
  % this speeds up the caching in real-time applications
  input.nSamples = 0;
end

current_argin = {input, desired};
if isequal(current_argin, previous_argin)
  % don't do the chantype detection again, but return the previous output from cache
  chantype = previous_argout{1};
  return
end

if isdata
  % the hdr, grad, elec or opto structure might have a different set of channels
  origlabel = input.label;

  if isfield(input, 'hdr')
    input = input.hdr;
    isheader = true;
  elseif isfield(input, 'grad')
    input = input.grad;
    isgrad = true;
  elseif isfield(input, 'elec')
    input = input.elec;
    iselec = true;
  elseif isfield(input, 'opto')
    input = input.opto;
    isopto = true;
  else
    % at least it contains channel labels
    islabel = true;
  end
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
  ft_error('the input that was provided to this function cannot be deciphered');
end

if isfield(input, 'chantype')
  % start with the provided channel types
  chantype = input.chantype(:);
else
  % start with unknown chantype for all channels
  chantype = repmat({'unknown'}, numchan, 1);
end

if ft_senstype(input, 'unknown')
  % don't bother doing all subsequent checks to determine the chantype of sensor array

elseif isheader && (ft_senstype(input, 'neuromag') || ft_senstype(input, 'babysquid74'))
  % channames-KI is the channel kind, 1=meg, 202=eog, 2=eeg, 3=trigger (I am not sure, but have inferred this from a single test file)
  % chaninfo-TY is the Coil chantype (0=magnetometer, 1=planar gradiometer)
  if isfield(input, 'orig') && isfield(input.orig, 'channames')
    for sel=find(input.orig.channames.KI(:)==202)'
      chantype{sel} = 'eog';
    end
    for sel=find(input.orig.channames.KI(:)==2)'
      chantype{sel} = 'eeg';
    end
    for sel=find(input.orig.channames.KI(:)==3)'
      chantype{sel} = 'digital trigger';
    end
    % determine the MEG channel subtype
    selmeg=find(input.orig.channames.KI(:)==1)';
    for i=1:length(selmeg)
      if input.orig.chaninfo.TY(i)==0
        chantype{selmeg(i)} = 'megmag';
      elseif input.orig.chaninfo.TY(i)==1
        % FIXME this might also be a axial gradiometer in case the BabySQUID data is read with the old reading routines
        chantype{selmeg(i)} = 'megplanar';
      end
    end

  elseif isfield(input, 'orig') && isfield(input.orig, 'chs') && isfield(input.orig.chs, 'coil_type')
    % all the chs.kinds and chs.coil_types are obtained from the MNE manual, p.210-211
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==2)' % planar gradiometers
      chantype(sel) = {'megplanar'}; %Neuromag-122 planar gradiometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3012)' %planar gradiometers
      chantype(sel) = {'megplanar'}; %Type T1 planar grad
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3013)' %planar gradiometers
      chantype(sel) = {'megplanar'}; %Type T2 planar grad
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3014)' %planar gradiometers
      chantype(sel) = {'megplanar'}; %Type T3 planar grad
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3022)' %magnetometers
      chantype(sel) = {'megmag'};    %Type T1 magenetometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3023)' %magnetometers
      chantype(sel) = {'megmag'};    %Type T2 magenetometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==3024)' %magnetometers
      chantype(sel) = {'megmag'};    %Type T3 magenetometer
    end
    for sel=find([input.orig.chs.kind]==1 & [input.orig.chs.coil_type]==7001)' %axial gradiometer
      chantype(sel) = {'megaxial'};
    end
    for sel=find([input.orig.chs.kind]==301)' %MEG reference channel, located far from head
      chantype(sel) = {'ref'};
    end
    for sel=find([input.orig.chs.kind]==2)'   %EEG channels
      chantype(sel) = {'eeg'};
    end
    for sel=find([input.orig.chs.kind]==201)' %MCG channels
      chantype(sel) = {'mcg'};
    end
    for sel=find([input.orig.chs.kind]==3)' %Stim channels
      if any([input.orig.chs(sel).logno] == 101) %new systems: 101 (and 102, if enabled) are digital; low numbers are 'pseudo-analog' (if enabled)
        chantype(sel([input.orig.chs(sel).logno] == 101)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] == 102)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] <= 32))  = {'analog trigger'};
        others = [input.orig.chs(sel).logno] > 32 & [input.orig.chs(sel).logno] ~= 101 & ...
          [input.orig.chs(sel).logno] ~= 102;
        chantype(sel(others)) = {'other trigger'};
      elseif any([input.orig.chs(sel).logno] == 14) %older systems: STI 014/015/016 are digital; lower numbers 'pseudo-analog'(if enabled)
        chantype(sel([input.orig.chs(sel).logno] == 14)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] == 15)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] == 16)) = {'digital trigger'};
        chantype(sel([input.orig.chs(sel).logno] <= 13)) = {'analog trigger'};
        others = [input.orig.chs(sel).logno] > 16;
        chantype(sel(others)) = {'other trigger'};
      else
        ft_warning('There does not seem to be a suitable trigger channel.');
        chantype(sel) = {'other trigger'};
      end
    end
    for sel=find([input.orig.chs.kind]==202)' %EOG
      chantype(sel) = {'eog'};
    end
    for sel=find([input.orig.chs.kind]==302)' %EMG
      chantype(sel) = {'emg'};
    end
    for sel=find([input.orig.chs.kind]==402)' %ECG
      chantype(sel) = {'ecg'};
    end
    for sel=find([input.orig.chs.kind]==502)' %MISC
      chantype(sel) = {'misc'};
    end
    for sel=find([input.orig.chs.kind]==602)' %Resp
      chantype(sel) = {'respiration'};
    end
  end

elseif ft_senstype(input, 'babysquid74')
  % the name can be something like "MEG 001" or "MEG001" or "MEG 0113" or "MEG0113"
  % i.e. with two or three digits and with or without a space
  sel = myregexp('^MEG', label);
  chantype(sel) = {'megaxial'};

elseif ft_senstype(input, 'neuromag122')
  % the name can be something like "MEG 001" or "MEG001" or "MEG 0113" or "MEG0113"
  % i.e. with two or three digits and with or without a space
  sel = myregexp('^MEG', label);
  chantype(sel) = {'megplanar'};

elseif ft_senstype(input, 'neuromag306') && isgrad
  % there should be 204 planar gradiometers and 102 axial magnetometers
  if isfield(input, 'tra')
    tmp = sum(abs(input.tra),2);
    sel = (tmp==median(tmp));
    chantype(sel) = {'megplanar'};
    sel = (tmp~=median(tmp));
    chantype(sel) = {'megmag'};
  end

elseif ft_senstype(input, 'neuromag306') && islabel
  sel = myregexp('^MEG.*1$', label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^MEG.*2$', label);
  chantype(sel) = {'megplanar'};
  sel = myregexp('^MEG.*3$', label);
  chantype(sel) = {'megplanar'};

elseif ft_senstype(input, 'neuromag306_combined') && islabel
  % the magnetometers are detected, the combined channels remain unknown
  sel = myregexp('^MEG.*1$', label);
  chantype(sel) = {'megmag'};

elseif ft_senstype(input, 'ctf') && isheader
  % The following is according to "CTF MEG(TM) File Formats" pdf, Release 5.2.1
  %
  % eMEGReference      0 Reference magnetometer channel
  % eMEGReference1     1 Reference 1st-order gradiometer channel
  % eMEGReference2     2 Reference 2nd-order gradiometer channel
  % eMEGReference3     3 Reference 3rd-order gradiometer channel
  % eMEGSensor         4 Sensor magnetometer channel located in head shell
  % eMEGSensor1        5 Sensor 1st-order gradiometer channel located in head shell
  % eMEGSensor2        6 Sensor 2nd-order gradiometer channel located in head shell
  % eMEGSensor3        7 Sensor 3rd-order gradiometer channel located in head shell
  % eEEGRef            8 EEG unipolar sensors not on the scalp
  % eEEGSensor         9 EEG unipolar sensors on the scalp
  % eADCRef           10 (see eADCAmpRef below)
  % eADCAmpRef        10 ADC amp channels from HLU or PIU (old electronics)
  % eStimRef          11 Stimulus channel for MEG41
  % eTimeRef          12 Time reference coming from video channel
  % ePositionRef      13 Not used
  % eDACRef           14 DAC channel from ECC or HLU
  % eSAMSensor        15 SAM channel derived through data analysis
  % eVirtualSensor    16 Virtual channel derived by combining two or more physical channels
  % eSystemTimeRef    17 System time showing elapsed time since trial started
  % eADCVoltRef       18 ADC volt channels from ECC
  % eStimAnalog       19 Analog trigger channels
  % eStimDigital      20 Digital trigger channels
  % eEEGBipolar       21 EEG bipolar sensor not on the scalp
  % eEEGAflg          22 EEG ADC over range flags
  % eMEGReset         23 MEG resets (counts sensor jumps for crosstalk purposes)
  % eDipSrc           24 Dipole source
  % eSAMSensorNorm    25 Normalized SAM channel derived through data analy- sis
  % eAngleRef         26 Orientation of head localization field
  % eExtractionRef    27 Extracted signal from each sensor of field generated by each localization coil
  % eFitErr           28 Fit error from each head localization coil
  % eOtherRef         29 Any other type of sensor not mentioned but still valid
  % eInvalidType      30 An invalid sensor

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
    ft_warning('could not determine channel chantype from the CTF header');
  end

  for sel=find(origSensType(:)==0)'
    chantype{sel} = 'refmag';
  end
  for sel=find(origSensType(:)==1)'
    chantype{sel} = 'refgrad';
  end
  for sel=find(origSensType(:)==5)'
    chantype{sel} = 'meggrad';
  end
  for sel=find(origSensType(:)==9)'
    chantype{sel} = 'eeg';
  end
  for sel=find(origSensType(:)==11)'
    % Stimulus channel for MEG41
    chantype{sel} = 'trigger';
  end
  for sel=find(origSensType(:)==13)'
    chantype{sel} = 'headloc'; % these represent the x, y, z position of the head coils
  end
  for sel=find(origSensType(:)==17)'
    chantype{sel} = 'clock';
  end
  for sel=find(origSensType(:)==18)'
    chantype{sel} = 'adc';
  end
  for sel=find(origSensType(:)==20)'
    % Digital trigger channels
    chantype{sel} = 'trigger';
  end
  for sel=find(origSensType(:)==28)'
    chantype{sel} = 'headloc_gof'; % these represent the goodness of fit for the head coils
  end
  for sel=find(origSensType(:)==29)'
    chantype{sel} = 'reserved'; % these are "reserved for future use", but relate to head localization
  end

elseif ft_senstype(input, 'ctf') && isgrad
  % in principle it is possible to look at the number of coils, but here the channels are identified based on their name
  sel = myregexp('^M[ZLR][A-Z][0-9][0-9]$', input.label);
  chantype(sel) = {'meggrad'};            % normal gradiometer channels
  sel = myregexp('^S[LR][0-9][0-9]$', input.label);
  chantype(sel) = {'meggrad'};            % normal gradiometer channels in the 64 channel CTF system
  sel = myregexp('^B[GPQR][0-9]$', input.label);
  chantype(sel) = {'refmag'};             % reference magnetometers
  sel = myregexp('^[GPQR][0-9][0-9]$', input.label);
  chantype(sel) = {'refgrad'};            % reference gradiometers

elseif ft_senstype(input, 'ctf') && islabel
  % the channels have to be identified based on their name alone
  sel = myregexp('^M[ZLR][A-Z][0-9][0-9]$', label);
  chantype(sel) = {'meggrad'};            % normal gradiometer channels
  sel = myregexp('^S[LR][0-9][0-9]$', label);
  chantype(sel) = {'meggrad'};            % normal gradiometer channels in the 64 channel CTF system
  sel = myregexp('^B[GPR][0-9]$', label);
  chantype(sel) = {'refmag'};             % reference magnetometers
  sel = myregexp('^[GPQR][0-9][0-9]$', label);
  chantype(sel) = {'refgrad'};            % reference gradiometers

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
    chantype(configtype==1 & numloops==1) = {'megmag'};
    chantype(configtype==1 & numloops==2) = {'meggrad'};
    chantype(configtype==2) = {'eeg'};
    chantype(configtype==3) = {'ref'}; % not known if mag or grad
    chantype(configtype==4) = {'aux'};
    chantype(configtype==5) = {'trigger'};

    % refine the distinction between refmag and refgrad to make the types
    % in grad and header consistent
    sel = myregexp('^M[CLR][xyz][aA]*$', label);
    chantype(sel) = {'refmag'};
    sel = myregexp('^G[xyz][xyz]A$', label);
    chantype(sel) = {'refgrad'};
  else
    % determine the chantype on the basis of the channel labels
    % all 4D-BTi MEG channels start with "A" followed by a number
    % all 4D-BTi reference channels start with M or G
    % all 4D-BTi EEG channels start with E, except for the 248-MEG/32-EEG system in Warsaw where they end with -1
    sel = myregexp('^A[0-9]+$', label);
    chantype(sel) = {'meg'};
    sel = myregexp('^M[CLR][xyz][aA]*$', label);
    chantype(sel) = {'refmag'};
    sel = myregexp('^G[xyz][xyz]A$', label);
    chantype(sel) = {'refgrad'};

    if isgrad && isfield(input, 'tra')
      gradtype = repmat({'unknown'}, size(input.label));
      gradtype(strncmp('A', input.label, 1)) = {'meg'};
      gradtype(strncmp('M', input.label, 1)) = {'refmag'};
      gradtype(strncmp('G', input.label, 1)) = {'refgrad'};
      % look at the number of coils of the meg channels
      selchan = find(strcmp('meg', gradtype));
      for k = 1:length(selchan)
        ncoils = length(find(input.tra(selchan(k),:)==1));
        if ncoils==1
          gradtype{selchan(k)} = 'megmag';
        elseif ncoils==2
          gradtype{selchan(k)} = 'meggrad';
        end
      end
      [selchan, selgrad] = match_str(label, input.label);
      chantype(selchan) = gradtype(selgrad);
    end

    % deal with additional channel types based on the names
    if isheader && issubfield(input, 'orig.channel_data.chan_label')
      tmplabel = {input.orig.channel_data.chan_label};
      tmplabel = tmplabel(:);
    else
      tmplabel = label; % might work
    end
    sel = find(strcmp('unknown', chantype));
    if ~isempty(sel)
      chantype(sel) = ft_chantype(tmplabel(sel));
      sel       = find(strcmp('unknown', chantype));
      if ~isempty(sel)
        % channels that start with E are assumed to be EEG
        % channels that end with -1 are also assumed to be EEG, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2389
        chantype(sel(cellfun(@(x) strcmp(x(end-1:end),'-1') || strcmp(x(1),'E'), label(sel)))) = {'eeg'};
      end
    end
  end

elseif ft_senstype(input, 'itab') && isheader
  origtype = [input.orig.ch.type];
  chantype(origtype==0) = {'unknown'};
  chantype(origtype==1) = {'ele'};
  chantype(origtype==2) = {'mag'}; % might be magnetometer or gradiometer, look at the number of coils
  chantype(origtype==4) = {'ele ref'};
  chantype(origtype==8) = {'mag ref'};
  chantype(origtype==16) = {'aux'};
  chantype(origtype==32) = {'param'};
  chantype(origtype==64) = {'digit'};
  chantype(origtype==128) = {'flag'};
  % these are the channels that are visible to FieldTrip
  chansel = 1:input.orig.nchan;
  chantype = chantype(chansel);

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
    label = input.label;
    sel = myregexp('[0-9][0-9][0-9]$', label);
    chantype(sel) = {'null'};
    sel = myregexp('^M[0-9][0-9][0-9]$', label);
    chantype(sel) = {'megmag'};
    sel = myregexp('^AG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'meggrad'};
    sel = myregexp('^PG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'megplanar'};
    sel = myregexp('^RM[0-9][0-9][0-9]$', label);
    chantype(sel) = {'refmag'};
    sel = myregexp('^RAG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'refgrad'};
    sel = myregexp('^RPG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'refplanar'};
    sel = myregexp('^TRIG[0-9][0-9][0-9]$', label);
    chantype(sel) = {'trigger'};
    %% Possible labels categorized in "eeg"
    sel_A = myregexp('^A[^G]*[0-9hzZ]$', label);
    sel_P = myregexp('^P[^G]*[0-9hzZ]$', label);
    sel_T = myregexp('^T[^R]*[0-9hzZ]$', label);
    sel_E = myregexp('^E$', label);
    sel_Z = myregexp('^[zZ]$', label);
    sel_M = myregexp('^M[0-9]$', label);
    sel_O = myregexp('^[BCFION]\w*[0-9hzZ]$', label);
    sel_EEG = myregexp('^EEG[0-9][0-9][0-9]$', label);
    sel = logical( sel_A + sel_P + sel_T + sel_E + sel_Z + sel_M + sel_O + sel_EEG );
    clear sel_A sel_P sel_T sel_E sel_Z sel_M sel_O sel_EEG
    chantype(sel) = {'eeg'};
    %% Additional EOG, ECG labels
    sel = myregexp('^EO[0-9]$', label); % EO
    chantype(sel) = {'eog'};
%    sel = myregexp('^ECG[0-9][0-9][0-9]$', label);
    sel_X = myregexp('^X[0-9]$', label); % X
    sel_ECG = myregexp('^ECG[0-9][0-9][0-9]$', label);
    sel = logical( sel_X + sel_ECG );
    clear sel_X sel_ECG
    chantype(sel) = {'ecg'};
    sel = myregexp('^ETC[0-9][0-9][0-9]$', label);
    chantype(sel) = {'etc'};

%   % shorten names
%    ch_info = input.orig.channel_info.channel;
%    type_orig = [ch_info.type];
%    sel = (type_orig == NullChannel);
%    chantype(sel) = {'null'};
%    sel = (type_orig == MagnetoMeter);
%    chantype(sel) = {'megmag'};
%    sel = (type_orig == AxialGradioMeter);
%    chantype(sel) = {'meggrad'};
%    sel = (type_orig == PlannerGradioMeter);
%    chantype(sel) = {'megplanar'};
%    sel = (type_orig == RefferenceMagnetoMeter);
%    chantype(sel) = {'refmag'};
%    sel = (type_orig == RefferenceAxialGradioMeter);
%    chantype(sel) = {'refgrad'};
%    sel = (type_orig == RefferencePlannerGradioMeter);
%    chantype(sel) = {'refplanar'};
%    sel = (type_orig == TriggerChannel);
%    chantype(sel) = {'trigger'};
%    sel = (type_orig == EegChannel);
%    chantype(sel) = {'eeg'};
%    sel = (type_orig == EcgChannel);
%    chantype(sel) = {'ecg'};
%    sel = (type_orig == EtcChannel);
%    chantype(sel) = {'etc'};

  elseif ft_hastoolbox('yokogawa')
    sel = (input.orig.channel_info(:, 2) == NullChannel);
    chantype(sel) = {'null'};
    sel = (input.orig.channel_info(:, 2) == MagnetoMeter);
    chantype(sel) = {'megmag'};
    sel = (input.orig.channel_info(:, 2) == AxialGradioMeter);
    chantype(sel) = {'meggrad'};
    sel = (input.orig.channel_info(:, 2) == PlannerGradioMeter);
    chantype(sel) = {'megplanar'};
    sel = (input.orig.channel_info(:, 2) == RefferenceMagnetoMeter);
    chantype(sel) = {'refmag'};
    sel = (input.orig.channel_info(:, 2) == RefferenceAxialGradioMeter);
    chantype(sel) = {'refgrad'};
    sel = (input.orig.channel_info(:, 2) == RefferencePlannerGradioMeter);
    chantype(sel) = {'refplanar'};
    sel = (input.orig.channel_info(:, 2) == TriggerChannel);
    chantype(sel) = {'trigger'};
    sel = (input.orig.channel_info(:, 2) == EegChannel);
    chantype(sel) = {'eeg'};
    sel = (input.orig.channel_info(:, 2) == EcgChannel);
    chantype(sel) = {'ecg'};
    sel = (input.orig.channel_info(:, 2) == EtcChannel);
    chantype(sel) = {'etc'};
  end

elseif ft_senstype(input, 'yokogawa') && isgrad
  % all channels in the gradiometer definition are meg
  % chantype(1:end) = {'meg'};
  % channels are identified based on their name: only magnetic as isgrad==1
  sel = myregexp('^M[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^AG[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'meggrad'};
  sel = myregexp('^PG[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'megplanar'};
  sel = myregexp('^RM[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'refmag'};
  sel = myregexp('^RAG[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'refgrad'};
  sel = myregexp('^RPG[0-9][0-9][0-9]$', input.label);
  chantype(sel) = {'refplanar'};

elseif ft_senstype(input, 'yokogawa') && islabel
  % the yokogawa channel labels are a mess, so autodetection is not possible
  % chantype(1:end) = {'meg'};
  sel = myregexp('[0-9][0-9][0-9]$', label);
  chantype(sel) = {'null'};
  sel = myregexp('^M[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^AG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'meggrad'};
  sel = myregexp('^PG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megplanar'};
  sel = myregexp('^RM[0-9][0-9][0-9]$', label);
  chantype(sel) = {'refmag'};
  sel = myregexp('^RAG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'refgrad'};
  sel = myregexp('^RPG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'refplanar'};
  sel = myregexp('^TRIG[0-9][0-9][0-9]$', label);
  chantype(sel) = {'trigger'};
  %% Possible labels categorized in "eeg"
  sel_A = myregexp('^A[^G]*[0-9hzZ]$', label);
  sel_P = myregexp('^P[^G]*[0-9hzZ]$', label);
  sel_T = myregexp('^T[^R]*[0-9hzZ]$', label);
  sel_E = myregexp('^E$', label);
  sel_Z = myregexp('^[zZ]$', label);
  sel_M = myregexp('^M[0-9]$', label);
  sel_O = myregexp('^[BCFION]\w*[0-9hzZ]$', label);
  sel_EEG = myregexp('^EEG[0-9][0-9][0-9]$', label);
  sel = logical( sel_A + sel_P + sel_T + sel_E + sel_Z + sel_M + sel_O + sel_EEG );
  clear sel_A sel_P sel_T sel_E sel_Z sel_M sel_O sel_EEG
  chantype(sel) = {'eeg'};
  %% Additional EOG, ECG labels
  sel = myregexp('^EO[0-9]$', label); % EO
  chantype(sel) = {'eog'};
% sel = myregexp('^ECG[0-9][0-9][0-9]$', label);
  sel_X = myregexp('^X[0-9]$', label); % X
  sel_ECG = myregexp('^ECG[0-9][0-9][0-9]$', label);
  sel = logical( sel_X + sel_ECG );
  clear sel_X sel_ECG
  chantype(sel) = {'ecg'};
  sel = myregexp('^ETC[0-9][0-9][0-9]$', label);
  chantype(sel) = {'etc'};

elseif ft_senstype(input, 'itab') && isheader
  sel = ([input.orig.ch.type]==0);
  chantype(sel) = {'unknown'};
  sel = ([input.orig.ch.type]==1);
  chantype(sel) = {'unknown'};
  sel = ([input.orig.ch.type]==2);
  chantype(sel) = {'megmag'};
  sel = ([input.orig.ch.type]==8);
  chantype(sel) = {'megref'};
  sel = ([input.orig.ch.type]==16);
  chantype(sel) = {'aux'};
  sel = ([input.orig.ch.type]==64);
  chantype(sel) = {'digital'};
  % not all channels are actually processed by FieldTrip, so only return
  % the types fopr the ones that read_header and read_data return
  chantype = chantype(input.orig.chansel);

elseif ft_senstype(input, 'itab') && isgrad
  % the channels have to be identified based on their name alone
  sel = myregexp('^MAG_[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^MAG_[0-9][0-9]$', label); % for the itab28 system
  chantype(sel) = {'megmag'};
  sel = myregexp('^MAG_[0-9]$', label); % for the itab28 system
  chantype(sel) = {'megmag'};
  sel = myregexp('^REF_[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megref'};
  sel = myregexp('^AUX.*$', label);
  chantype(sel) = {'aux'};

elseif ft_senstype(input, 'itab') && islabel
  % the channels have to be identified based on their name alone
  sel = myregexp('^MAG_[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megmag'};
  sel = myregexp('^REF_[0-9][0-9][0-9]$', label);
  chantype(sel) = {'megref'};
  sel = myregexp('^AUX.*$', label);
  chantype(sel) = {'aux'};

elseif ft_senstype(input, 'eeg') && islabel
  % use an external helper function to define the list with EEG channel names
  chantype(match_str(label, ft_senslabel('eeg1005'))) = {'eeg'};          % this includes all channels from the 1010 and 1020 arrangement
  chantype(match_str(label, ft_senslabel(ft_senstype(input)))) = {'eeg'}; % this will work for biosemi, egi and other detected channel arrangements

elseif ft_senstype(input, 'eeg') && iselec
  % all channels in an electrode definition must be eeg channels
  chantype(:) = {'eeg'};

elseif ft_senstype(input, 'eeg') && isheader
  % use an external helper function to define the list with EEG channel names
  chantype(match_str(input.label, ft_senslabel(ft_senstype(input)))) = {'eeg'};

elseif ft_senstype(input, 'plexon') && isheader
  % this is a complete header that was read from a Plexon *.nex file using read_plexon_nex
  for i=1:numchan
    switch input.orig.VarHeader(i).Type
      case 0
        chantype{i} = 'spike';
      case 1
        chantype{i} = 'event';
      case 2
        chantype{i} = 'interval';  % Interval variables?
      case 3
        chantype{i} = 'waveform';
      case 4
        chantype{i} = 'population'; % Population variables ?
      case 5
        chantype{i} = 'analog';
      otherwise
        % keep the default 'unknown' chantype
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
    chantype(intersect(strmatch(label2type{i}{j}, lower(label)), find(strcmp(chantype, 'unknown')))) = label2type{i}(1);
  end
end

if isdata
  % the input was replaced by one of hdr, grad, elec, opto
  [sel1, sel2] = match_str(origlabel, input.label);
  origtype = repmat({'unknown'}, size(sel1));
  origtype(sel1) = chantype(sel2);
  % the hdr, grad, elec or opto structure might have a different set of channels
  chantype = origtype;
end

if all(strcmp(chantype, 'unknown')) && ~recursion
  % try whether only lowercase channel labels makes a difference
  if islabel
    recursion = true;
    chantype = ft_chantype(lower(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = lower(input.label);
    recursion = true;
    chantype = ft_chantype(input);
    recursion = false;
  end
end

if all(strcmp(chantype, 'unknown')) && ~recursion
  % try whether only uppercase channel labels makes a difference
  if islabel
    recursion = true;
    chantype = ft_chantype(upper(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = upper(input.label);
    recursion = true;
    chantype = ft_chantype(input);
    recursion = false;
  end
end

if nargin>1
  % return a boolean vector
  if isequal(desired, 'meg') || isequal(desired, 'ref')
    % only compare the first three characters, i.e. meggrad or megmag should match
    chantype = strncmp(desired, chantype, 3);
  else
    chantype = strcmp(desired, chantype);
  end
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {chantype};
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
