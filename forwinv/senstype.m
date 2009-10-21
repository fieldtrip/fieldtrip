function [type] = senstype(input, desired)

% SENSTYPE determines the type of sensors by looking at the channel names
% and comparing them with predefined lists.
%
% Use as
%   [type] = senstype(sens)
% to get a string describing the type, or
%   [flag] = senstype(sens, desired)
% to get a boolean value.
%
% The output type can be any of the following
%   'electrode'
%   'magnetometer'
%   'biosemi64'
%   'biosemi128'
%   'biosemi256'
%   'bti148'
%   'bti148_planar'
%   'bti248'
%   'bti248_planar'
%   'ctf151'
%   'ctf151_planar'
%   'ctf275'
%   'ctf275_planar'
%   'egi128'
%   'egi256'
%   'egi32'
%   'egi64'
%   'ext1020'
%   'neuromag122'
%   'neuromag306'
%   'yokogawa160'
%   'yokogawa160_planar'
%   'plexon'
%   'itab153'
%   'itab153_planar'
%
% The optional input argument for the desired type can be any of the above,
% or any of the following
%   'eeg'
%   'meg'
%   'meg_planar'
%   'meg_axial'
%   'ctf'
%   'bti'
%   'neuromag'
%   'yokogawa'
%
% Besides specifiying a grad or elec structure as input, also allowed is
% giving a data structure containing a grad or elec field, or giving a list
% of channel names (as cell-arrray). I.e. assuming a FieldTrip data
% structure, all of the following calls would be correct.
%   senstype(data)
%   senstype(data.label)
%   senstype(data.grad)
%   senstype(data.grad.label)
%
% See also SENSLABEL, CHANTYPE, READ_SENS, COMPUTE_LEADFIELD 

% Copyright (C) 2007-2009, Robert Oostenveld
%
% $Log: senstype.m,v $
% Revision 1.23  2009/10/19 13:13:32  roboos
% cleaned up autodetection of the type of input, it should return 'unknown' on unexpected inputs (e.g. headshapes) and not give an error
%
% Revision 1.22  2009/10/19 10:12:52  vlalit
% Fixed a bug that caused shape to be recognized as elec.
%
% Revision 1.21  2009/10/16 12:27:53  roboos
% some small changes pertaining to the itab/chieti format
%
% Revision 1.20  2009/10/13 10:39:19  roboos
% added support for 153 channel itab system
%
% Revision 1.19  2009/07/29 08:04:38  roboos
% cleaned up the code, no functional change
%
% Revision 1.18  2009/07/28 11:16:23  roboos
% removed keyboard statement, thanks to Jurrian
%
% Revision 1.17  2009/07/28 10:17:41  roboos
% make distinction between only label, label+pnt and label+pnt+ori
%
% Revision 1.16  2009/07/27 16:04:51  roboos
% improved distinction between eeg and meg, fixes problem with biosemi-eeg being detected as "ctf" due to reference channel match
%
% Revision 1.15  2009/06/19 16:51:50  vlalit
% Added biosemi64 system of  Diane Whitmer, I don't know how generic it is.
%
% Revision 1.14  2009/05/07 13:34:09  roboos
% added ctf64
%
% Revision 1.13  2009/04/01 06:51:43  roboos
% implemented caching for the type detection in case the same input is given multiple times
% this uses a persistent variable
%
% Revision 1.12  2009/03/06 08:50:23  roboos
% added handling of plexon header
%
% Revision 1.11  2009/02/02 16:27:41  roboos
% changed order of detecting sens.grad/elec on Vladimirs request, don't know why
%
% Revision 1.10  2008/09/10 09:12:11  roboos
% added alternative definition of channel names without a space in the label for neuromag 122 and 306
%
% Revision 1.9  2008/09/10 07:53:27  roboos
% moved definition of channel label sets to seperate function
%
% Revision 1.8  2008/03/18 13:40:27  roboos
% added quick fix for ctf275_planar, the problem is that one channel is missing from the ctf275 list (i.e. it is only 274 long)
%
% Revision 1.7  2008/03/18 12:25:06  roboos
% use sens.type if available, this requires that the content of sens.type is consistent with the strings returned by this function
% preallocate cell-arrays
%
% Revision 1.6  2008/02/29 15:25:31  roboos
% fixed bug for eeg (thanks to Doug), updated documentation
%
% Revision 1.5  2008/02/29 14:04:42  roboos
% added gradiometers to btiref
% added general types to the checks at the end
%
% Revision 1.4  2008/02/29 13:52:12  roboos
% added bti248 and planar
% changed order of arguments for ismember, needed when a lot of EEG channels are present in ctf151
% changed order in which checks are performed, first ctf275 and only then ctf151
%
% Revision 1.3  2008/02/28 09:16:35  roboos
% fixed bug due to many trigger and headloc channels in ctf275 labels
%
% Revision 1.2  2008/02/27 17:01:54  roboos
% added a whole list of channel names, fixed some bugs
%
% Revision 1.1  2007/07/25 08:31:12  roboos
% implemented new helper function
%

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

if nargin<2
  % ensure that all input arguments are defined
  desired = [];
end

current_argin = {input, desired};
if isequal(current_argin, previous_argin)
  % don't do the type detection again, but return the previous values from cache
  type = previous_argout{1};
  return
end

isdata   = isa(input, 'struct') && isfield(input, 'hdr')   && isfield(input.hdr, 'label');
isheader = isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'Fs');
isgrad   = isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'pnt')  &&  isfield(input, 'ori');
iselec   = isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'pnt')  && ~isfield(input, 'ori');
islabel  = isa(input, 'cell')   && isa(input{1}, 'char');

% the input may be a data structure which then contains a grad/elec structure
if isdata
  if isfield(input.hdr, 'grad')
    sens = input.hdr.grad;
    isgrad = true;
  elseif isfield(input.hdr, 'elec')
    sens = input.hdr.elec;
    iselec = true;
  else
    sens.label = input.hdr.label;
    islabel = true;
  end
elseif isheader
  if isfield(input, 'grad')
    sens = input.grad;
    isgrad = true;
  elseif isfield(input, 'elec')
    sens   = input.elec;
    iselec = true;
  else
    sens.label = input.label;
    islabel = true;
  end
elseif isgrad
  sens = input;
elseif iselec
  sens = input;
elseif islabel
  sens.label = input;
else
  sens = [];
end

if isfield(sens, 'type')
  % preferably the structure specifies its own type
  type = sens.type;

elseif issubfield(sens, 'orig.FileHeader') &&  issubfield(sens, 'orig.VarHeader')
  % this is a complete header that was read from a Plexon *.nex file using read_plexon_nex
  type = 'plexon';

else
  % start with unknown, then try to determine the proper type by looking at the labels
  type = 'unknown';

  if isgrad
    % probably this is MEG, determine the type of magnetometer/gradiometer system
    % note that the order here is important: first check whether it matches a 275 channel system, then a 151 channel system, since the 151 channels are a subset of the 275
    if     (mean(ismember(senslabel('ctf275'),        sens.label)) > 0.8)
      type = 'ctf275';
    elseif (mean(ismember(senslabel('ctfheadloc'),    sens.label)) > 0.8)  % look at the head localization channels
      type = 'ctf275';
    elseif (mean(ismember(senslabel('ctf151'),        sens.label)) > 0.8)
      type = 'ctf151';
    elseif (mean(ismember(senslabel('ctf64'),         sens.label)) > 0.8)
      type = 'ctf64';
    elseif (mean(ismember(senslabel('ctf275_planar'), sens.label)) > 0.8)
      type = 'ctf275_planar';
    elseif (mean(ismember(senslabel('ctf151_planar'), sens.label)) > 0.8)
      type = 'ctf151_planar';
    elseif (mean(ismember(senslabel('bti248'),        sens.label)) > 0.8)
      type = 'bti248';
    elseif (mean(ismember(senslabel('bti148'),        sens.label)) > 0.8)
      type = 'bti148';
    elseif (mean(ismember(senslabel('bti248_planar'), sens.label)) > 0.8)
      type = 'bti248_planar';
    elseif (mean(ismember(senslabel('bti148_planar'), sens.label)) > 0.8)
      type = 'bti148_planar';
    elseif (mean(ismember(senslabel('itab153'),       sens.label)) > 0.8)
      type = 'itab153';
    elseif (mean(ismember(senslabel('itab153_planar'), sens.label)) > 0.8)
      type = 'itab153_planar';
    elseif (mean(ismember(senslabel('neuromag306'),   sens.label)) > 0.8)
      type = 'neuromag306';
    elseif (mean(ismember(senslabel('neuromag306alt'),sens.label)) > 0.8)  % an alternative set without spaces in the name
      type = 'neuromag306';
    elseif (mean(ismember(senslabel('neuromag122'),   sens.label)) > 0.8)
      type = 'neuromag122';
    elseif (mean(ismember(senslabel('neuromag122alt'),sens.label)) > 0.8)  % an alternative set without spaces in the name
      type = 'neuromag122';
    elseif any(ismember(senslabel('btiref'), sens.label))
      type = 'bti'; % it might be 148 or 248 channels
    elseif any(ismember(senslabel('ctfref'), sens.label))
      type = 'ctf'; % it might be 151 or 275 channels
    elseif isfield(sens, 'pnt') && isfield(sens, 'ori') && numel(sens.label)==numel(sens.pnt)
      type = 'magnetometer';
    else
      type = 'meg';
    end

  elseif iselec
    % probably this is EEG
    if     (mean(ismember(senslabel('biosemi256'),    sens.label)) > 0.8)
      type = 'biosemi256';
    elseif (mean(ismember(senslabel('biosemi128'),    sens.label)) > 0.8)
      type = 'biosemi128';
    elseif (mean(ismember(senslabel('biosemi64'),     sens.label)) > 0.8)
      type = 'biosemi64';
    elseif (mean(ismember(senslabel('egi256'),        sens.label)) > 0.8)
      type = 'egi256';
    elseif (mean(ismember(senslabel('egi128'),        sens.label)) > 0.8)
      type = 'egi128';
    elseif (mean(ismember(senslabel('egi64'),         sens.label)) > 0.8)
      type = 'egi64';
    elseif (mean(ismember(senslabel('egi32'),         sens.label)) > 0.8)
      type = 'egi32';
    elseif (sum(ismember(sens.label,         senslabel('eeg1005'))) > 10) % Otherwise it's not even worth recognizing
      type = 'ext1020';
    else
      type = 'electrode';
    end

  elseif islabel
    % look only at the channel labels
    if     (mean(ismember(senslabel('ctf275'),        sens.label)) > 0.8)
      type = 'ctf275';
    elseif (mean(ismember(senslabel('ctfheadloc'),    sens.label)) > 0.8)  % look at the head localization channels
      type = 'ctf275';
    elseif (mean(ismember(senslabel('ctf151'),        sens.label)) > 0.8)
      type = 'ctf151';
    elseif (mean(ismember(senslabel('ctf64'),         sens.label)) > 0.8)
      type = 'ctf64';
    elseif (mean(ismember(senslabel('ctf275_planar'), sens.label)) > 0.8)
      type = 'ctf275_planar';
    elseif (mean(ismember(senslabel('ctf151_planar'), sens.label)) > 0.8)
      type = 'ctf151_planar';
    elseif (mean(ismember(senslabel('bti248'),        sens.label)) > 0.8)
      type = 'bti248';
    elseif (mean(ismember(senslabel('bti148'),        sens.label)) > 0.8)
      type = 'bti148';
    elseif (mean(ismember(senslabel('bti248_planar'), sens.label)) > 0.8)
      type = 'bti248_planar';
    elseif (mean(ismember(senslabel('bti148_planar'), sens.label)) > 0.8)
      type = 'bti148_planar';
    elseif (mean(ismember(senslabel('itab153'),       sens.label)) > 0.8)
      type = 'itab153';
    elseif (mean(ismember(senslabel('itab153_planar'), sens.label)) > 0.8)
      type = 'itab153_planar';
    elseif (mean(ismember(senslabel('neuromag306'),   sens.label)) > 0.8)
      type = 'neuromag306';
    elseif (mean(ismember(senslabel('neuromag306alt'),sens.label)) > 0.8)  % an alternative set without spaces in the name
      type = 'neuromag306';
    elseif (mean(ismember(senslabel('neuromag122'),   sens.label)) > 0.8)
      type = 'neuromag122';
    elseif (mean(ismember(senslabel('neuromag122alt'),sens.label)) > 0.8)  % an alternative set without spaces in the name
      type = 'neuromag122';
    elseif (mean(ismember(senslabel('biosemi256'),    sens.label)) > 0.8)
      type = 'biosemi256';
    elseif (mean(ismember(senslabel('biosemi128'),    sens.label)) > 0.8)
      type = 'biosemi128';
    elseif (mean(ismember(senslabel('biosemi64'),     sens.label)) > 0.8)
      type = 'biosemi64';
    elseif (mean(ismember(senslabel('egi256'),        sens.label)) > 0.8)
      type = 'egi256';
    elseif (mean(ismember(senslabel('egi128'),        sens.label)) > 0.8)
      type = 'egi128';
    elseif (mean(ismember(senslabel('egi64'),         sens.label)) > 0.8)
      type = 'egi64';
    elseif (mean(ismember(senslabel('egi32'),         sens.label)) > 0.8)
      type = 'egi32';
    elseif (sum(ismember(sens.label,         senslabel('eeg1005'))) > 10) % Otherwise it's not even worth recognizing
      type = 'ext1020';
    elseif any(ismember(senslabel('btiref'), sens.label))
      type = 'bti'; % it might be 148 or 248 channels
    elseif any(ismember(senslabel('ctfref'), sens.label))
      type = 'ctf'; % it might be 151 or 275 channels
    end

  end % look at label, ori and/or pnt
end % if isfield(sens, 'type')

if ~isempty(desired)
  % return a boolean flag
  switch desired
    case 'eeg'
      type = any(strcmp(type, {'eeg' 'electrode' 'biosemi64' 'biosemi128' 'biosemi256' 'egi32' 'egi64' 'egi128' 'egi256' 'ext1020'}));
    case 'biosemi'
      type = any(strcmp(type, {'biosemi64' 'biosemi128' 'biosemi256'}));
    case 'egi'
      type = any(strcmp(type, {'egi64' 'egi128' 'egi256'}));
    case 'meg'
      type = any(strcmp(type, {'meg' 'magnetometer' 'ctf' 'bti' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar' 'neuromag122' 'neuromag306' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar' 'yokogawa160' 'yokogawa160_planar'}));
    case 'ctf'
      type = any(strcmp(type, {'ctf' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar'}));
    case 'bti'
      type = any(strcmp(type, {'bti' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar'}));
    case 'neuromag'
      type = any(strcmp(type, {'neuromag122' 'neuromag306'}));
    case 'yokogawa'
      type = any(strcmp(type, {'yokogawa160' 'yokogawa160_planar'}));
    case 'itab'
      type = any(strcmp(type, {'itab153' 'itab153_planar'}));
    case 'meg_axial'
      % note that neuromag306 is mixed planar and axial
      type = any(strcmp(type, {'magnetometer' 'neuromag306' 'ctf151' 'ctf275' 'bti148' 'bti248' 'yokogawa160'}));
    case 'meg_planar'
      % note that neuromag306 is mixed planar and axial
      type = any(strcmp(type, {'neuromag122' 'neuromag306' 'ctf151_planar' 'ctf275_planar' 'bti148_planar' 'bti248_planar' 'yokogawa160_planar'}));
    otherwise
      type = any(strcmp(type, desired));
  end % switch desired
end % detemine the correspondence to the desired type

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {type};
if isempty(previous_argin)
  previous_argin  = current_argin;
  previous_argout = current_argout;
end

return % senstype main()
