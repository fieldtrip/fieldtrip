function [type] = ft_senstype(input, desired)

% FT_SENSTYPE determines the type of sensors by looking at the channel names
% and comparing them with predefined lists.
%
% Use as
%   [type] = ft_senstype(sens)
% to get a string describing the type, or
%   [flag] = ft_senstype(sens, desired)
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
%   ft_senstype(data)
%   ft_senstype(data.label)
%   ft_senstype(data.grad)
%   ft_senstype(data.grad.label)
%
% See also FT_SENSLABEL, FT_CHANTYPE, FT_READ_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2007-2009, Robert Oostenveld
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
% $Id: ft_senstype.m 2075 2010-11-05 12:03:49Z roboos $

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

if iscell(input) && numel(input)<4 && ~all(cellfun(@ischar, input))
  % this might represent combined EEG, ECoG and/or MEG
  type = cell(size(input));
  if nargin<2
    desired = cell(size(input)); % empty elements
  end
  for i=1:numel(input)
    type{i} = ft_senstype(input{i}, desired{i});
  end
  return
end

if nargin<2
  % ensure that all input arguments are defined
  desired = [];
end

current_argin = {input, desired};
if isequal(current_argin, previous_argin)
  % don't do the type detection again, but return the previous values from
  % cache
  type = previous_argout{1};
  return
end

% FIXME the detection of the type of input structure should perhaps be done using the datatype function
isdata   = isa(input, 'struct') && isfield(input, 'hdr');
isheader = isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'Fs');
isgrad   = isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'pnt')  &&  isfield(input, 'ori');
iselec   = isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'pnt')  && ~isfield(input, 'ori');
islabel  = isa(input, 'cell')   && ~isempty(input) && isa(input{1}, 'char');

if ~isdata && ~isheader
  % timelock or freq structures don't have the header structure
  % the header is also removed from raw data after megrealign
  % the gradiometer definition is lost after megplanar+combineplanar
  isdata = isa(input, 'struct') && (isfield(input, 'grad') || isfield(input, 'elec') || isfield(input, 'label'));
end

% the input may be a data structure which then contains a grad/elec structure, a header or only the labels
if isdata
  % preferably look at the data and not the header for the grad, because it might be re-balanced and/or planar
  if isfield(input, 'grad')
    sens = input.grad;
    isgrad = true;
  elseif issubfield(input, 'hdr.grad')
    sens = input.hdr.grad;
    isgrad = true;
  elseif issubfield(input, 'hdr.elec')
    sens = input.hdr.elec;
    iselec = true;
  elseif isfield(input, 'elec')
    sens = input.elec;
    iselec = true;
  elseif issubfield(input, 'hdr.label')
    sens.label = input.hdr.label;
    islabel = true;
  elseif isfield(input, 'label')
    sens.label = input.label;
    islabel = true;
  end
elseif isheader
  if isfield(input, 'grad')
    sens = input.grad;
    isgrad = true;
  elseif isfield(input, 'elec')
    sens   = input.elec;
    iselec = true;
  elseif isfield(input, 'label')
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

if isfield(input, 'type')
  % preferably the structure specifies its own type
  type = input.type;

elseif issubfield(input, 'orig.FileHeader') &&  issubfield(input, 'orig.VarHeader')
  % this is a complete header that was read from a Plexon *.nex file using read_plexon_nex
  type = 'plexon';

elseif issubfield(input, 'orig.stname')
  % this is a complete header that was read from an ITAB dataset
  type = 'itab';

elseif issubfield(input, 'orig.sys_name')
  % this is a complete header that was read from a Yokogawa dataset
  type = 'yokogawa160';

elseif issubfield(input, 'orig.FILE.Ext') && strcmp(input.orig.FILE.Ext, 'edf')
  % this is a complete header that was read from an EDF or EDF+ dataset
  type = 'electrode';

else
  % start with unknown, then try to determine the proper type by looking at the labels
  type = 'unknown';

  if isgrad
    % probably this is MEG, determine the type of magnetometer/gradiometer system
    % note that the order here is important: first check whether it matches a 275 channel system, then a 151 channel system, since the 151 channels are a subset of the 275
    if     (mean(ismember(ft_senslabel('ctf275'),        sens.label)) > 0.8)
      type = 'ctf275';
    elseif (mean(ismember(ft_senslabel('ctfheadloc'),    sens.label)) > 0.8)  % look at the head localization channels
      type = 'ctf275';
    elseif (mean(ismember(ft_senslabel('ctf151'),        sens.label)) > 0.8)
      type = 'ctf151';
    elseif (mean(ismember(ft_senslabel('ctf64'),         sens.label)) > 0.8)
      type = 'ctf64';
    elseif (mean(ismember(ft_senslabel('ctf275_planar'), sens.label)) > 0.8)
      type = 'ctf275_planar';
    elseif (mean(ismember(ft_senslabel('ctf151_planar'), sens.label)) > 0.8)
      type = 'ctf151_planar';
    elseif (mean(ismember(ft_senslabel('bti248'),        sens.label)) > 0.8)
      type = 'bti248';
    elseif (mean(ismember(ft_senslabel('bti148'),        sens.label)) > 0.8)
      type = 'bti148';
    elseif (mean(ismember(ft_senslabel('bti248_planar'), sens.label)) > 0.8)
      type = 'bti248_planar';
    elseif (mean(ismember(ft_senslabel('bti148_planar'), sens.label)) > 0.8)
      type = 'bti148_planar';
    elseif (mean(ismember(ft_senslabel('itab153'),       sens.label)) > 0.8)
      type = 'itab153';
    elseif (mean(ismember(ft_senslabel('itab153_planar'), sens.label)) > 0.8)
      type = 'itab153_planar';
    elseif (mean(ismember(ft_senslabel('yokogawa160'),    sens.label)) > 0.4)
      type = 'yokogawa160';
    elseif (mean(ismember(ft_senslabel('yokogawa160_planar'), sens.label)) > 0.4)
      type = 'yokogawa160_planar';
    elseif (mean(ismember(ft_senslabel('neuromag306'),   sens.label)) > 0.8)
      type = 'neuromag306';
    elseif (mean(ismember(ft_senslabel('neuromag306alt'),sens.label)) > 0.8)  % an alternative set without spaces in the name
      type = 'neuromag306';
    elseif (mean(ismember(ft_senslabel('neuromag122'),   sens.label)) > 0.8)
      type = 'neuromag122';
    elseif (mean(ismember(ft_senslabel('neuromag122alt'),sens.label)) > 0.8)  % an alternative set without spaces in the name
      type = 'neuromag122';
    elseif any(ismember(ft_senslabel('btiref'), sens.label))
      type = 'bti'; % it might be 148 or 248 channels
    elseif any(ismember(ft_senslabel('ctfref'), sens.label))
      type = 'ctf'; % it might be 151 or 275 channels
    elseif isfield(sens, 'pnt') && isfield(sens, 'ori') && numel(sens.label)==size(sens.pnt,1)
      warning('could be Yokogawa system');
      type = 'magnetometer';
    else
      warning('could be Yokogawa system');
      type = 'meg';
    end

  elseif iselec
    % probably this is EEG
    if     (mean(ismember(ft_senslabel('biosemi256'),    sens.label)) > 0.8)
      type = 'biosemi256';
    elseif (mean(ismember(ft_senslabel('biosemi128'),    sens.label)) > 0.8)
      type = 'biosemi128';
    elseif (mean(ismember(ft_senslabel('biosemi64'),     sens.label)) > 0.8)
      type = 'biosemi64';
    elseif (mean(ismember(ft_senslabel('egi256'),        sens.label)) > 0.8)
      type = 'egi256';
    elseif (mean(ismember(ft_senslabel('egi128'),        sens.label)) > 0.8)
      type = 'egi128';
    elseif (mean(ismember(ft_senslabel('egi64'),         sens.label)) > 0.8)
      type = 'egi64';
    elseif (mean(ismember(ft_senslabel('egi32'),         sens.label)) > 0.8)
      type = 'egi32';
    elseif (sum(ismember(sens.label,         ft_senslabel('eeg1005'))) > 10) % Otherwise it's not even worth recognizing
      type = 'ext1020';
    else
      type = 'electrode';
    end

  elseif islabel
    % look only at the channel labels
    if     (mean(ismember(ft_senslabel('ctf275'),        sens.label)) > 0.8)
      type = 'ctf275';
    elseif (mean(ismember(ft_senslabel('ctfheadloc'),    sens.label)) > 0.8)  % look at the head localization channels
      type = 'ctf275';
    elseif (mean(ismember(ft_senslabel('ctf151'),        sens.label)) > 0.8)
      type = 'ctf151';
    elseif (mean(ismember(ft_senslabel('ctf64'),         sens.label)) > 0.8)
      type = 'ctf64';
    elseif (mean(ismember(ft_senslabel('ctf275_planar'), sens.label)) > 0.8)
      type = 'ctf275_planar';
    elseif (mean(ismember(ft_senslabel('ctf151_planar'), sens.label)) > 0.8)
      type = 'ctf151_planar';
    elseif (mean(ismember(ft_senslabel('bti248'),        sens.label)) > 0.8)
      type = 'bti248';
    elseif (mean(ismember(ft_senslabel('bti148'),        sens.label)) > 0.8)
      type = 'bti148';
    elseif (mean(ismember(ft_senslabel('bti248_planar'), sens.label)) > 0.8)
      type = 'bti248_planar';
    elseif (mean(ismember(ft_senslabel('bti148_planar'), sens.label)) > 0.8)
      type = 'bti148_planar';
    elseif (mean(ismember(ft_senslabel('itab153'),       sens.label)) > 0.8)
      type = 'itab153';
    elseif (mean(ismember(ft_senslabel('itab153_planar'), sens.label)) > 0.8)
      type = 'itab153_planar';
    elseif (mean(ismember(ft_senslabel('yokogawa160'),    sens.label)) > 0.4)
      type = 'yokogawa160';
    elseif (mean(ismember(ft_senslabel('yokogawa160_planar'), sens.label)) > 0.4)
      type = 'yokogawa160_planar';
    elseif any(mean(ismember(ft_senslabel('neuromag306'),   sens.label)) > 0.8)
      type = 'neuromag306';
    elseif any(mean(ismember(ft_senslabel('neuromag306alt'),sens.label)) > 0.8)  % an alternative set without spaces in the name
      type = 'neuromag306';
    elseif any(mean(ismember(ft_senslabel('neuromag122'),   sens.label)) > 0.8)
      type = 'neuromag122';
    elseif any(mean(ismember(ft_senslabel('neuromag122alt'),sens.label)) > 0.8)  % an alternative set without spaces in the name
      type = 'neuromag122';
    elseif (mean(ismember(ft_senslabel('biosemi256'),    sens.label)) > 0.8)
      type = 'biosemi256';
    elseif (mean(ismember(ft_senslabel('biosemi128'),    sens.label)) > 0.8)
      type = 'biosemi128';
    elseif (mean(ismember(ft_senslabel('biosemi64'),     sens.label)) > 0.8)
      type = 'biosemi64';
    elseif (mean(ismember(ft_senslabel('egi256'),        sens.label)) > 0.8)
      type = 'egi256';
    elseif (mean(ismember(ft_senslabel('egi128'),        sens.label)) > 0.8)
      type = 'egi128';
    elseif (mean(ismember(ft_senslabel('egi64'),         sens.label)) > 0.8)
      type = 'egi64';
    elseif (mean(ismember(ft_senslabel('egi32'),         sens.label)) > 0.8)
      type = 'egi32';
    elseif (sum(ismember(sens.label,         ft_senslabel('eeg1005'))) > 10) % Otherwise it's not even worth recognizing
      type = 'ext1020';
    elseif any(ismember(ft_senslabel('btiref'), sens.label))
      type = 'bti'; % it might be 148 or 248 channels
    elseif any(ismember(ft_senslabel('ctfref'), sens.label))
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
      type = any(strcmp(type, {'itab' 'itab153' 'itab153_planar'}));
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

return % ft_senstype main()
