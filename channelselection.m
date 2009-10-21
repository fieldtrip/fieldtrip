function [channel] = channelselection(channel, datachannel)

% CHANNELSELECTION for EEG and MEG labels
%
% This function translates the user-specified list of channels into channel
% labels as they occur in the data. This channel selection procedure can
% be used throughout fieldtrip.
%
% You can specify a mixture of real channel labels and of special strings,
% or index numbers that will be replaced by the corresponding channel
% labels. Channels that are not present in the raw datafile are
% automatically removed from the channel list.
%
% E.g.
%  'all'     is replaced by all channels in the datafile
%  'gui'     a graphical user interface will pop up to select the channels
%  'C*'      is replaced by all channels that match the wildcard, e.g. C1, C2, C3, ...
%  '*1'      is replaced by all channels that match the wildcard, e.g. C1, P1, F1, ...
%  'M*1'     is replaced by all channels that match the wildcard, e.g. MEG0111, MEG0131, MEG0131, ...
%  'MEG'     is replaced by all MEG channels (works for CTF, 4D and Neuromag)
%  'MEGREF'  is replaced by all MEG reference channels (works for CTF and 4D)
%  'EEG'     is replaced by all recognized EEG channels (this is system dependent)
%  'EEG1020' is replaced by 'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', ...
%  'EOG'     is replaced by all recognized EOG channels
%  'EMG'     is replaced by all channels in the datafile starting with 'EMG'
%  'lfp'     is replaced by all channels in the datafile starting with
%  'lfp'
%  'mua'     is replaced by all channels in the datafile starting with 'mua'
%  'spike'   is replaced by all channels in the datafile starting with 'spike'
%  10        is replaced by the 10th channel in the datafile
%
% Other channel groups are
%   'EEG1010'    with approximately 90 electrodes
%   'EEG1005'    with approximately 350 electrodes
%   'EEGCHWILLA' for Dorothee Chwilla's electrode caps (used at the NICI)
%   'EEGBHAM'    for the 128 channel EEG system used in Birmingham
%   'EEGREF'     for mastoid and ear electrodes (M1, M2, LM, RM, A1, A2)
%   'MZ'         for MEG zenith
%   'ML'         for MEG left
%   'MR'         for MEG right
%   'MLx', 'MRx' and 'MZx' with x=C,F,O,P,T for left/right central, frontal, occipital, parietal and temporal
%
% You can also exclude channels or channel groups using the following syntax
%   {'all', '-POz', '-Fp1', -EOG'}

% Note that the order of channels that is returned should correspond with
% the order of the channels in the data.

% Copyright (C) 2003-2009, Robert Oostenveld
%
% $Log: channelselection.m,v $
% Revision 1.39  2009/08/04 13:54:20  roboos
% allow datachannel as single string
%
% Revision 1.38  2009/07/08 07:27:37  roboos
% improved wildcard selection, now also in middle like "M*1"
%
% Revision 1.37  2009/06/19 16:51:30  vlalit
% Added biosemi64 system of  Diane Whitmer, I don't know how generic it is.
%
% Revision 1.36  2009/02/18 07:57:07  roboos
% added support for selection based on wildcards (like 'C*' and '*1')
%
% Revision 1.35  2009/01/21 16:59:08  jansch
% added possibility to select (sets of) references for 4D data
%
% Revision 1.34  2008/11/10 15:16:22  roboos
% added snippet of code to speed up the channel selection if there is a perfect match, no need to look for channel groups and immediately bail out
%
% Revision 1.33  2008/10/13 12:29:07  jansch
% added plain 'bti' to known 4D senstypes
%
% Revision 1.32  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.31  2008/09/10 09:11:35  roboos
% added definition of EEG for neuromag
%
% Revision 1.30  2008/09/10 08:39:33  roboos
% make stronger use of senstype and senslabel helper functions
% now also supports egi and biosemi systems for EEG (thanks to Vladimir)
%
% Revision 1.29  2008/08/13 14:03:52  roboos
% added hEOG and vEOG
%
% Revision 1.28  2008/07/17 14:45:39  roboos
% document megref, added megref for ctf
%
% Revision 1.27  2008/07/17 10:25:20  roboos
% ensure that the output ordering is the same as in teh original list with channel names
%
% Revision 1.26  2008/05/21 09:47:14  jansch
% added MEGREF (for 4d-systems) as a channel-group
%
% Revision 1.25  2008/04/21 14:38:26  jansch
% added support 'MEG' for NM122 systems
%
% Revision 1.24  2008/02/27 17:03:02  roboos
% use senstype to distinguish between ctf and bti
%
% Revision 1.23  2008/01/31 20:12:03  roboos
% It is possible given the naming structure of the channels in the
% 4D file that both the isctf and isbti flags can be set in the script
% on lines 155 and 156.  An extra IF statement has been added to
% detect and correct this situation in favour of the isbti flag.
% [thanks to Gavin]
%

fieldtripdefs

if any(size(channel) == 0)
  % there is nothing to do if it is empty
  return
end

if isnumeric(channel)
  % change index into channelname
  channel = datachannel(channel);
  return
end

if ~iscell(channel)
  % ensure that a single input argument like 'all' also works
  channel = {channel};
end

if ~iscell(datachannel)
  % ensure that a single input argument like 'all' also works
  datachannel = {datachannel};
end

% ensure that both inputs are column vectors
channel     = channel(:);
datachannel = datachannel(:);

% remove channels that occur more than once, this sorts the channels alphabetically
[channel, indx] = unique(channel);
% undo the sorting, make the order identical to that of the data channels
[dum, indx] = sort(indx);
channel = channel(indx);

[dataindx, chanindx] = match_str(datachannel, channel);
if length(chanindx)==length(channel)
  % there is a perfect match between the channels and the datachannels, only some reordering is needed
  channel = channel(chanindx);
  % no need to look at channel groups
  return
end

% define the known groups with channel labels
labelall  = datachannel;
label1020 = senslabel('eeg1020'); % use external helper function
label1010 = senslabel('eeg1010'); % use external helper function
label1005 = senslabel('eeg1005'); % use external helper function
labelchwilla = {'Fz', 'Cz', 'Pz', 'F7', 'F8', 'LAT', 'RAT', 'LT', 'RT', 'LTP', 'RTP', 'OL', 'OR', 'FzA', 'Oz', 'F7A', 'F8A', 'F3A', 'F4A', 'F3', 'F4', 'P3', 'P4', 'T5', 'T6', 'P3P', 'P4P'}';
labelbham = {'P9', 'PPO9h', 'PO7', 'PPO5h', 'PPO3h', 'PO5h', 'POO9h', 'PO9', 'I1', 'OI1h', 'O1', 'POO1', 'PO3h', 'PPO1h', 'PPO2h', 'POz', 'Oz', 'Iz', 'I2', 'OI2h', 'O2', 'POO2', 'PO4h', 'PPO4h', 'PO6h', 'POO10h', 'PO10', 'PO8', 'PPO6h', 'PPO10h', 'P10', 'P8', 'TPP9h', 'TP7', 'TTP7h', 'CP5', 'TPP7h', 'P7', 'P5', 'CPP5h', 'CCP5h', 'CP3', 'P3', 'CPP3h', 'CCP3h', 'CP1', 'P1', 'Pz', 'CPP1h', 'CPz', 'CPP2h', 'P2', 'CPP4h', 'CP2', 'CCP4h', 'CP4', 'P4', 'P6', 'CPP6h', 'CCP6h', 'CP6', 'TPP8h', 'TP8', 'TPP10h', 'T7', 'FTT7h', 'FT7', 'FC5', 'FCC5h', 'C5', 'C3', 'FCC3h', 'FC3', 'FC1', 'C1', 'CCP1h', 'Cz', 'FCC1h', 'FCz', 'FFC1h', 'Fz', 'FFC2h', 'FC2', 'FCC2h', 'CCP2h', 'C2', 'C4', 'FCC4h', 'FC4', 'FC6', 'FCC6h', 'C6', 'TTP8h', 'T8', 'FTT8h', 'FT8', 'FT9', 'FFT9h', 'F7', 'FFT7h', 'FFC5h', 'F5', 'AFF7h', 'AF7', 'AF5h', 'AFF5h', 'F3', 'FFC3h', 'F1', 'AF3h', 'Fp1', 'Fpz', 'Fp2', 'AFz', 'AF4h', 'F2', 'FFC4h', 'F4', 'AFF6h', 'AF6h', 'AF8', 'AFF8h', 'F6', 'FFC6h', 'FFT8h', 'F8', 'FFT10h', 'FT10'};
labelref  = {'M1', 'M2', 'LM', 'RM', 'A1', 'A2'}';
labeleog  = datachannel(strncmp('EOG', datachannel, length('EOG')));               % anything that starts with EOG
labeleog  = {labeleog{:} 'HEOG', 'VEOG', 'VEOG-L', 'VEOG-R', 'hEOG', 'vEOG'}';     % or any of these
labelemg  = datachannel(strncmp('EMG', datachannel, length('EMG')));
labellfp  = datachannel(strncmp('lfp', datachannel, length('lfp')));
labelmua  = datachannel(strncmp('mua', datachannel, length('mua')));
labelspike  = datachannel(strncmp('spike', datachannel, length('spike')));

% use regular expressions to deal with the wildcards
labelreg = false(size(datachannel));
findreg = [];
for i=1:length(channel)
  if length(channel{i})>1 && channel{i}(1)=='*'
    % the wildcard is at the start
    labelreg = labelreg | ~cellfun(@isempty, regexp(datachannel, ['.*' channel{i}(2:end) '$'], 'once'));
    findreg  = [findreg i];
  elseif length(channel{i})>1 && channel{i}(end)=='*'
    % the wildcard is at the end
    labelreg = labelreg | ~cellfun(@isempty, regexp(datachannel, ['^' channel{i}(1:end-1) '.*'], 'once'));
    findreg  = [findreg i];
  elseif length(channel{i})>1 && any(channel{i}=='*')
    % the wildcard is in the middle
    sel  = strfind(channel{i}, '*');
    str1 = channel{i}(1:(sel-1));
    str2 = channel{i}((sel+1):end);
    labelreg = labelreg | ~cellfun(@isempty, regexp(datachannel, ['^' str1 '.*' str2 '$'], 'once'));
    findreg  = [findreg i];
  end
end
labelreg = datachannel(labelreg);

% initialize all the system-specific variables to empty
labelmeg   = [];
labelmref  = [];
labeleeg   = [];

switch senstype(datachannel)

  case {'ctf', 'ctf275', 'ctf151', 'ctf275_planar', 'ctf151_planar'}
    % all CTF MEG channels start with "M"
    % all CTF reference channels start with B, G, P, Q or R
    % all CTF EEG channels start with "EEG"
    labelmeg = datachannel(strncmp('M'  , datachannel, length('M'  )));
    labelmref = [datachannel(strncmp('B'  , datachannel, 1));
      datachannel(strncmp('G'  , datachannel, 1));
      datachannel(strncmp('P'  , datachannel, 1));
      datachannel(strncmp('Q'  , datachannel, 1));
      datachannel(strncmp('R'  , datachannel, length('G'  )))];
    labeleeg  = datachannel(strncmp('EEG', datachannel, length('EEG')));

    % Not sure whether this should be here or outside the switch or
    % whether these specifications should be supported for systems
    % other than CTF.
    labelmz   = datachannel(strncmp('MZ' , datachannel, length('MZ' )));	% central MEG channels
    labelml   = datachannel(strncmp('ML' , datachannel, length('ML' )));	% left    MEG channels
    labelmr   = datachannel(strncmp('MR' , datachannel, length('MR' )));	% right   MEG channels
    labelmlc  = datachannel(strncmp('MLC', datachannel, length('MLC')));
    labelmlf  = datachannel(strncmp('MLF', datachannel, length('MLF')));
    labelmlo  = datachannel(strncmp('MLO', datachannel, length('MLO')));
    labelmlp  = datachannel(strncmp('MLP', datachannel, length('MLP')));
    labelmlt  = datachannel(strncmp('MLT', datachannel, length('MLT')));
    labelmrc  = datachannel(strncmp('MRC', datachannel, length('MRC')));
    labelmrf  = datachannel(strncmp('MRF', datachannel, length('MRF')));
    labelmro  = datachannel(strncmp('MRO', datachannel, length('MRO')));
    labelmrp  = datachannel(strncmp('MRP', datachannel, length('MRP')));
    labelmrt  = datachannel(strncmp('MRT', datachannel, length('MRT')));
    labelmzc  = datachannel(strncmp('MZC', datachannel, length('MZC')));
    labelmzf  = datachannel(strncmp('MZF', datachannel, length('MZF')));
    labelmzo  = datachannel(strncmp('MZO', datachannel, length('MZO')));
    labelmzp  = datachannel(strncmp('MZP', datachannel, length('MZP')));

  case {'bti', 'bti248', 'bti148', 'bti248_planar', 'bti148_planar'}
    % all 4D-BTi MEG channels start with "A"
    % all 4D-BTi reference channels start with M or G
    labelmeg = datachannel(strncmp('A'  , datachannel, 1));
    labelmref = [datachannel(strncmp('M'  , datachannel, 1));
      datachannel(strncmp('G'  , datachannel, 1))];
    labelmrefa = datachannel(~cellfun(@isempty,strfind(datachannel, 'a')));
    labelmrefc = datachannel(strncmp('MC', datachannel, 2));
    labelmrefg = datachannel(strncmp('G', datachannel,  1));
    labelmrefl = datachannel(strncmp('ML', datachannel, 2));
    labelmrefr = datachannel(strncmp('MR', datachannel, 2));

  case {'neuromag306', 'neuromag122'}
    % all neuromag MEG channels start with MEG
    % all neuromag EEG channels start with EEG
    labelmeg = datachannel(strncmp('MEG', datachannel, length('MEG')));
    labeleeg = datachannel(strncmp('EEG', datachannel, length('EEG')));

  case {'biosemi64', 'biosemi128', 'biosemi256', 'egi64', 'egi128', 'egi256', 'ext1020'}
    % use an external helper function to define the list with EEG channel names
    labeleeg = senslabel(senstype(datachannel));

end % switch senstype

% figure out if there are bad channels or channel groups that should be excluded
findbadchannel = strncmp('-', channel, length('-'));      % bad channels start with '-'
badchannel = channel(findbadchannel);
if ~isempty(badchannel)
  for i=1:length(badchannel)
    badchannel{i} = badchannel{i}(2:end);                 % remove the '-' from the channel label
  end
  badchannel = channelselection(badchannel, datachannel); % support exclusion of channel groups
  channel(findbadchannel) = [];                           % remove them from the channels to be processed
end

% determine if any of the known groups is mentioned in the channel list
findall        = find(strcmp(channel, 'all'));
% findreg (for the wildcards) is dealt with in the channel group specification above
findmeg        = find(strcmp(channel, 'MEG'));
findemg        = find(strcmp(channel, 'EMG'));
findeeg        = find(strcmp(channel, 'EEG'));
findeeg1020    = find(strcmp(channel, 'EEG1020'));
findeeg1010    = find(strcmp(channel, 'EEG1010'));
findeeg1005    = find(strcmp(channel, 'EEG1005'));
findeegchwilla = find(strcmp(channel, 'EEGCHWILLA'));
findeegbham    = find(strcmp(channel, 'EEGBHAM'));
findeegref     = find(strcmp(channel, 'EEGREF'));
findmegref     = find(strcmp(channel, 'MEGREF'));
findmegrefa    = find(strcmp(channel, 'MEGREFA'));
findmegrefc    = find(strcmp(channel, 'MEGREFC'));
findmegrefg    = find(strcmp(channel, 'MEGREFG'));
findmegrefl    = find(strcmp(channel, 'MEGREFL'));
findmegrefr    = find(strcmp(channel, 'MEGREFR'));
findeog        = find(strcmp(channel, 'EOG'));
findmz         = find(strcmp(channel, 'MZ' ));
findml         = find(strcmp(channel, 'ML' ));
findmr         = find(strcmp(channel, 'MR' ));
findmlc        = find(strcmp(channel, 'MLC'));
findmlf        = find(strcmp(channel, 'MLF'));
findmlo        = find(strcmp(channel, 'MLO'));
findmlp        = find(strcmp(channel, 'MLP'));
findmlt        = find(strcmp(channel, 'MLT'));
findmrc        = find(strcmp(channel, 'MRC'));
findmrf        = find(strcmp(channel, 'MRF'));
findmro        = find(strcmp(channel, 'MRO'));
findmrp        = find(strcmp(channel, 'MRP'));
findmrt        = find(strcmp(channel, 'MRT'));
findmzc        = find(strcmp(channel, 'MZC'));
findmzf        = find(strcmp(channel, 'MZF'));
findmzo        = find(strcmp(channel, 'MZO'));
findmzp        = find(strcmp(channel, 'MZP'));
findlfp        = find(strcmp(channel, 'lfp'));
findmua        = find(strcmp(channel, 'mua'));
findspike      = find(strcmp(channel, 'spike'));
findgui        = find(strcmp(channel, 'gui'));

% remove any occurences of groups in the channel list
channel([
  findall
  findreg
  findmeg
  findemg
  findeeg
  findeeg1020
  findeeg1010
  findeeg1005
  findeegchwilla
  findeegbham
  findeegref
  findmegref
  findeog
  findmz
  findml
  findmr
  findmlc
  findmlf
  findmlo
  findmlp
  findmlt
  findmrc
  findmrf
  findmro
  findmrp
  findmrt
  findmzc
  findmzf
  findmzo
  findmzp
  findlfp
  findmua
  findspike
  findgui
  ]) = [];

% add the full channel labels to the channel list
if findall,        channel = [channel; labelall]; end
if findreg,        channel = [channel; labelreg]; end
if findmeg,        channel = [channel; labelmeg]; end
if findemg,        channel = [channel; labelemg]; end
if findeeg,        channel = [channel; labeleeg]; end
if findeeg1020,    channel = [channel; label1020]; end
if findeeg1010,    channel = [channel; label1010]; end
if findeeg1005,    channel = [channel; label1005]; end
if findeegchwilla, channel = [channel; labelchwilla]; end
if findeegbham,    channel = [channel; labelbham]; end
if findeegref,     channel = [channel; labelref]; end
if findmegref,     channel = [channel; labelmref]; end
if findmegrefa,    channel = [channel; labelmrefa]; end
if findmegrefc,    channel = [channel; labelmrefc]; end
if findmegrefg,    channel = [channel; labelmrefg]; end
if findmegrefl,    channel = [channel; labelmrefl]; end
if findmegrefr,    channel = [channel; labelmrefr]; end
if findeog,        channel = [channel; labeleog]; end
if findmz ,        channel = [channel; labelmz ]; end
if findml ,        channel = [channel; labelml ]; end
if findmr ,        channel = [channel; labelmr ]; end
if findmlc,        channel = [channel; labelmlc]; end
if findmlf,        channel = [channel; labelmlf]; end
if findmlo,        channel = [channel; labelmlo]; end
if findmlp,        channel = [channel; labelmlp]; end
if findmlt,        channel = [channel; labelmlt]; end
if findmrc,        channel = [channel; labelmrc]; end
if findmrf,        channel = [channel; labelmrf]; end
if findmro,        channel = [channel; labelmro]; end
if findmrp,        channel = [channel; labelmrp]; end
if findmrt,        channel = [channel; labelmrt]; end
if findmzc,        channel = [channel; labelmzc]; end
if findmzf,        channel = [channel; labelmzf]; end
if findmzo,        channel = [channel; labelmzo]; end
if findmzp,        channel = [channel; labelmzp]; end
if findlfp,        channel = [channel; labellfp]; end
if findmua,        channel = [channel; labelmua]; end
if findspike,      channel = [channel; labelspike]; end

% remove channel labels that have been excluded by the user
badindx = match_str(channel, badchannel);
channel(badindx) = [];

% remove channel labels that are not present in the data
chanindx = match_str(channel, datachannel);

channel  = channel(chanindx);

if findgui
  indx = select_channel_list(datachannel, match_str(datachannel, channel), 'Select channels');
  channel = datachannel(indx);
end

% remove channels that occur more than once, this sorts the channels alphabetically
channel = unique(channel);

% undo the sorting, make the order identical to that of the data channels
[dataindx, indx] = match_str(datachannel, channel);
channel = channel(indx);
