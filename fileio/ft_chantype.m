function type = ft_chantype(input, desired)

% FT_CHANTYPE determines for each channel what type it is, e.g. planar/axial gradiometer or magnetometer
%
% Use as
%   type = ft_chantype(hdr)
%   type = ft_chantype(sens)
%   type = ft_chantype(label)
% or as
%   type = ft_chantype(hdr, desired)
%   type = ft_chantype(sens, desired)
%   type = ft_chantype(label, desired)

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% determine the type of input

isheader = isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'Fs');
isgrad   = isa(input, 'struct') && isfield(input, 'pnt') && isfield(input, 'ori');
islabel  = isa(input, 'cell')   && isa(input{1}, 'char');

hdr   = input;
grad  = input;
label = input;

if isheader
  numchan = length(hdr.label);
  if isfield(hdr, 'grad')
    grad  = hdr.grad;
  end
  label = hdr.label;
elseif isgrad
  label   = grad.label;
  numchan = length(label);
elseif islabel
  numchan = length(label);
else
  error('the input that was provided to this function cannot be deciphered');
end

% start with unknown type
type = cell(numchan,1);
for i=1:length(type)
  type{i} = 'unknown';
end

if senstype(input, 'neuromag')
  % channames-KI is the channel kind, 1=meg, 202=eog, 2=eeg, 3=trigger (I am nut sure, but have inferred this from a single test file)
  % chaninfo-TY is the Coil type (0=magnetometer, 1=planar gradiometer)
  if isfield(hdr, 'orig') && isfield(hdr.orig, 'channames')
    for sel=find(hdr.orig.channames.KI(:)==202)'
      type{sel} = 'eog';
    end
    for sel=find(hdr.orig.channames.KI(:)==2)'
      type{sel} = 'eeg';
    end
    for sel=find(hdr.orig.channames.KI(:)==3)'
      type{sel} = 'digital trigger';
    end
    % determinge the MEG channel subtype
    selmeg=find(hdr.orig.channames.KI(:)==1)';
    for i=1:length(selmeg)
      if hdr.orig.chaninfo.TY(i)==0
        type{selmeg(i)} = 'megmag';
      elseif hdr.orig.chaninfo.TY(i)==1
        type{selmeg(i)} = 'megplanar';
      end
    end

  elseif isfield(hdr, 'orig') && isfield(hdr.orig, 'chs') && isfield(hdr.orig.chs, 'coil_type')
    %all the chs.kinds and chs.coil_types are obtained from the MNE
    %manual, p.210-211
    for sel=find([hdr.orig.chs.kind]==1 & [hdr.orig.chs.coil_type]==2)' %planar gradiometers
      type(sel) = {'megplanar'}; %Neuromag-122 planar gradiometer
    end
    for sel=find([hdr.orig.chs.kind]==1 & [hdr.orig.chs.coil_type]==3012)' %planar gradiometers
      type(sel) = {'megplanar'}; %Type T1 planar grad
    end
    for sel=find([hdr.orig.chs.kind]==1 & [hdr.orig.chs.coil_type]==3013)' %planar gradiometers
      type(sel) = {'megplanar'}; %Type T2 planar grad
    end
    for sel=find([hdr.orig.chs.kind]==1 & [hdr.orig.chs.coil_type]==3014)' %planar gradiometers
      type(sel) = {'megplanar'}; %Type T3 planar grad
    end
    for sel=find([hdr.orig.chs.kind]==1 & [hdr.orig.chs.coil_type]==3022)' %magnetometers
      type(sel) = {'megmag'};    %Type T1 magenetometer
    end
    for sel=find([hdr.orig.chs.kind]==1 & [hdr.orig.chs.coil_type]==3023)' %magnetometers
      type(sel) = {'megmag'};    %Type T2 magenetometer
    end
    for sel=find([hdr.orig.chs.kind]==1 & [hdr.orig.chs.coil_type]==3024)' %magnetometers
      type(sel) = {'megmag'};    %Type T3 magenetometer
    end
    for sel=find([hdr.orig.chs.kind]==301)' %MEG reference channel, located far from head
      type(sel) = {'ref'};
    end
    for sel=find([hdr.orig.chs.kind]==2)' %EEG channels
      type(sel) = {'eeg'};
    end
    for sel=find([hdr.orig.chs.kind]==201)' %MCG channels
      type(sel) = {'mcg'};
    end
    for sel=find([hdr.orig.chs.kind]==3)' %Stim channels
      if any([hdr.orig.chs(sel).logno] == 101) %new systems: 101 (and 102, if enabled) are digital;
        %low numbers are 'pseudo-analog' (if enabled)
        type(sel([hdr.orig.chs(sel).logno] == 101)) = {'digital trigger'};
        type(sel([hdr.orig.chs(sel).logno] == 102)) = {'digital trigger'};
        type(sel([hdr.orig.chs(sel).logno] <= 32)) = {'analog trigger'};
        others = [hdr.orig.chs(sel).logno] > 32 & [hdr.orig.chs(sel).logno] ~= 101 & ...
          [hdr.orig.chs(sel).logno] ~= 102;
        type(sel(others)) = {'other trigger'};
      elseif any([hdr.orig.chs(sel).logno] == 14) %older systems: STI 014/015/016 are digital;
        %lower numbers 'pseudo-analog'(if enabled)
        type(sel([hdr.orig.chs(sel).logno] == 14)) = {'digital trigger'};
        type(sel([hdr.orig.chs(sel).logno] == 15)) = {'digital trigger'};
        type(sel([hdr.orig.chs(sel).logno] == 16)) = {'digital trigger'};
        type(sel([hdr.orig.chs(sel).logno] <= 13)) = {'analog trigger'};
        others = [hdr.orig.chs(sel).logno] > 16;
        type(sel(others)) = {'other trigger'};
      else
        warning('There does not seem to be a suitable trigger channel.');
        type(sel) = {'other trigger'};
      end
    end
    for sel=find([hdr.orig.chs.kind]==202)' %EOG
      type(sel) = {'eog'};
    end
    for sel=find([hdr.orig.chs.kind]==302)' %EMG
      type(sel) = {'emg'};
    end
    for sel=find([hdr.orig.chs.kind]==402)' %ECG
      type(sel) = {'ecg'};
    end
    for sel=find([hdr.orig.chs.kind]==502)'  %MISC
      type(sel) = {'misc'};
    end
    for sel=find([hdr.orig.chs.kind]==602)' %Resp
      type(sel) = {'respiration'};
    end
  end

elseif senstype(input, 'ctf') && isheader
  % meg channels are 5, refmag 0, refgrad 1, adcs 18, trigger 11, eeg 9
  if isfield(hdr, 'orig') && isfield(hdr.orig, 'sensType')
    origSensType = hdr.orig.sensType;
  elseif isfield(hdr, 'orig') && isfield(hdr.orig, 'res4')
    origSensType =  [hdr.orig.res4.senres.sensorTypeIndex];
  else
    warning('could not determine channel type from the CTF header');
    origSensType = [];
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

elseif senstype(input, 'ctf') && isgrad
  % in principle it is possible to look at the number of coils, but here the channels are identified based on their name
  sel = myregexp('^M[ZLR][A-Z][0-9][0-9]$', grad.label);
  type(sel) = {'meggrad'};            % normal gradiometer channels
  sel = myregexp('^B[GPR][0-9]$', grad.label);
  type(sel) = {'refmag'};             % reference magnetometers
  sel = myregexp('^[GPQR][0-9][0-9]$', grad.label);
  type(sel) = {'refgrad'};            % reference gradiometers

elseif senstype(input, 'ctf') && islabel
  % the channels have to be identified based on their name alone
  sel = myregexp('^M[ZLR][A-Z][0-9][0-9]$', label);
  type(sel) = {'meggrad'};                % normal gradiometer channels
  sel = myregexp('^B[GPR][0-9]$', label);
  type(sel) = {'refmag'};             % reference magnetometers
  sel = myregexp('^[GPQR][0-9][0-9]$', label);
  type(sel) = {'refgrad'};            % reference gradiometers

elseif senstype(input, 'bti')
  % all 4D-BTi MEG channels start with "A"
  % all 4D-BTi reference channels start with M or G
  % all 4D-BTi EEG channels start with E
  type(strncmp('A', label, 1)) = {'meg'};
  type(strncmp('M', label, 1)) = {'refmag'};
  type(strncmp('G', label, 1)) = {'refgrad'};

  if isfield(grad, 'tra')
    selchan = find(strcmp('mag', type));
    for k = 1:length(selchan)
      ncoils = length(find(grad.tra(selchan(k),:)==1));
      if ncoils==1,
        type{selchan(k)} = 'megmag';
      elseif ncoils==2,
        type{selchan(k)} = 'meggrad';
      end
    end
  end

  % This is to allow setting additional channel types based on the names
  if isheader && issubfield(hdr, 'orig.channel_data.chan_label')
    tmplabel = {hdr.orig.channel_data.chan_label};
    tmplabel = tmplabel(:);
  else
    tmplabel = label; % might work
  end
  sel      = strmatch('unknown', type, 'exact');
  if ~isempty(sel)
    type(sel)= ft_chantype(tmplabel(sel));
    sel      = strmatch('unknown', type, 'exact');
    if ~isempty(sel)
      type(sel(strncmp('E', label(sel), 1))) = {'eeg'};
    end
  end

elseif senstype(input, 'yokogawa') && isheader
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
    sel = (hdr.orig.channel_info(:, 2) == NullChannel);
    type(sel) = {'null'};
    sel = (hdr.orig.channel_info(:, 2) == MagnetoMeter);
    type(sel) = {'megmag'};
    sel = (hdr.orig.channel_info(:, 2) == AxialGradioMeter);
    type(sel) = {'meggrad'};
    sel = (hdr.orig.channel_info(:, 2) == PlannerGradioMeter);
    type(sel) = {'megplanar'};
    sel = (hdr.orig.channel_info(:, 2) == RefferenceMagnetoMeter);
    type(sel) = {'refmag'};
    sel = (hdr.orig.channel_info(:, 2) == RefferenceAxialGradioMeter);
    type(sel) = {'refgrad'};
    sel = (hdr.orig.channel_info(:, 2) == RefferencePlannerGradioMeter);
    type(sel) = {'refplanar'};
    sel = (hdr.orig.channel_info(:, 2) == TriggerChannel);
    type(sel) = {'trigger'};
    sel = (hdr.orig.channel_info(:, 2) == EegChannel);
    type(sel) = {'eeg'};
    sel = (hdr.orig.channel_info(:, 2) == EcgChannel);
    type(sel) = {'ecg'};
    sel = (hdr.orig.channel_info(:, 2) == EtcChannel);
    type(sel) = {'etc'};

elseif senstype(input, 'yokogawa') && isgrad
  % all channels in the gradiometer definition are meg
  type(1:end) = {'meg'};
  
elseif senstype(input, 'yokogawa') && islabel
  % the yokogawa channel labels are a mess, so autodetection is not possible
  type(1:end) = {'meg'};

elseif senstype(input, 'itab') && isheader
  sel = ([hdr.orig.ch.type]==0);
  type(sel) = {'unknown'};
  sel = ([hdr.orig.ch.type]==1);
  type(sel) = {'unknown'};
  sel = ([hdr.orig.ch.type]==2);
  type(sel) = {'megmag'};
  sel = ([hdr.orig.ch.type]==8);
  type(sel) = {'megref'};
  sel = ([hdr.orig.ch.type]==16);
  type(sel) = {'aux'};
  sel = ([hdr.orig.ch.type]==64);
  type(sel) = {'digital'};
  % not all channels are actually processed by fieldtrip, so only return
  % the types fopr the ones that read_header and read_data return
  type = type(hdr.orig.chansel);

elseif senstype(input, 'itab') && isgrad
  % the channels have to be identified based on their name alone
  sel = myregexp('^MAG_[0-9][0-9][0-9]$', label);
  type(sel) = {'megmag'};
  sel = myregexp('^REF_[0-9][0-9][0-9]$', label);
  type(sel) = {'megref'};
  sel = myregexp('^AUX.*$', label);
  type(sel) = {'aux'};

elseif senstype(input, 'itab') && islabel
  % the channels have to be identified based on their name alone
  sel = myregexp('^MAG_[0-9][0-9][0-9]$', label);
  type(sel) = {'megmag'};
  sel = myregexp('^REF_[0-9][0-9][0-9]$', label);
  type(sel) = {'megref'};
  sel = myregexp('^AUX.*$', label);
  type(sel) = {'aux'};

elseif senstype(input, 'eeg') && islabel
  % use an external helper function to define the list with EEG channel names
  type(match_str(label, senslabel(senstype(label)))) = {'eeg'};

elseif senstype(input, 'eeg') && isheader
  % use an external helper function to define the list with EEG channel names
  type(match_str(hdr.label, senslabel(senstype(hdr)))) = {'eeg'};

elseif senstype(input, 'plexon') && isheader
  % this is a complete header that was read from a Plexon *.nex file using read_plexon_nex
  for i=1:numchan
    switch hdr.orig.VarHeader(i).Type
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

end % senstype

% if possible, set additional types based on channel labels
label2type = {
  {'ecg', 'ekg'};
  {'emg'};
  {'eog', 'heog', 'veog'};
  {'lfp'};
  {'eeg'};
  };
for i = 1:numel(label2type)
  for j = 1:numel(label2type{i})
    type(intersect(strmatch(label2type{i}{j}, lower(label)),...
      strmatch('unknown', type, 'exact'))) = label2type{i}(1);
  end
end

if nargin>1
  % return a boolean vector
  if isequal(desired, 'meg') || isequal(desired, 'ref')
    type = strncmp(desired, type, 3);
  else
    type = strcmp(desired, type);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function match = myregexp(pat, list)
match = false(size(list));
for i=1:numel(list)
  match(i) = ~isempty(regexp(list{i}, pat, 'once'));
end


