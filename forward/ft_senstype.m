function [type] = ft_senstype(input, desired)

% FT_SENSTYPE determines the type of acquisition device by looking at the
% channel names and comparing them with predefined lists.
%
% Use as
%   [type] = ft_senstype(sens)
%   [flag] = ft_senstype(sens, desired)
%
% The output type can be any of the following
%   'ctf64'
%   'ctf151'
%   'ctf151_planar'
%   'ctf275'
%   'ctf275_planar'
%   'bti148'
%   'bti148_planar'
%   'bti248'
%   'bti248_planar'
%   'bti248grad'
%   'bti248grad_planar'
%   'itab28'
%   'itab153'
%   'itab153_planar'
%   'yokogawa9'
%   'yokogawa64'
%   'yokogawa64_planar'
%   'yokogawa160'
%   'yokogawa160_planar'
%   'yokogawa440'
%   'ext1020'       this includes eeg1020, eeg1010 and eeg1005)
%   'neuromag122'
%   'neuromag306'
%   'babysquid74'   this is a BabySQUID system from Tristan Technologies
%   'artemis123'    this is a BabySQUID system from Tristan Technologies
%   'magview'       this is a BabySQUID system from Tristan Technologies
%   'egi32'
%   'egi64'
%   'egi128'
%   'egi256'
%   'biosemi64'
%   'biosemi128'
%   'biosemi256'
%   'neuralynx'
%   'plexon'
%   'artinis'
%   'eeg' (this was called 'electrode' in older versions)
%   'meg' (this was called 'magnetometer' in older versions)
%   'nirs'
%
% The optional input argument for the desired type can be any of the above,
% or any of the following generic classes of acquisition systems
%   'eeg'
%   'meg'
%   'meg_planar'
%   'meg_axial'
%   'ctf'
%   'bti'
%   'neuromag'
%   'yokogawa'
%   'babysquid'
% If you specify the desired type, this function will return a boolean
% true/false depending on the input data.
%
% Besides specifiying a sensor definition (i.e. a grad or elec structure,
% see FT_DATATYPE_SENS), it is also possible to give a data structure
% containing a grad or elec field, or giving a list of channel names (as
% cell-arrray). So assuming that you have a FieldTrip data structure, any
% of the following calls would also be fine.
%   ft_senstype(hdr)
%   ft_senstype(data)
%   ft_senstype(data.label)
%   ft_senstype(data.grad)
%   ft_senstype(data.grad.label)
%
% See also FT_SENSLABEL, FT_CHANTYPE, FT_READ_SENS, FT_COMPUTE_LEADFIELD, FT_DATATYPE_SENS

% Copyright (C) 2007-2015, Robert Oostenveld
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
% $Id$

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

% this is to avoid a recursion loop
persistent recursion
if isempty(recursion)
  recursion = false;
end

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
  % don't do the type detection again, but return the previous output from cache
  type = previous_argout{1};
  return
end

isdata   = isa(input, 'struct')  && (isfield(input, 'hdr') || isfield(input, 'time') || isfield(input, 'freq') || isfield(input, 'grad') || isfield(input, 'elec') || isfield(input, 'opto'));
isheader = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'Fs');
isgrad   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  &&  isfield(input, 'ori'); % old style
iselec   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  && ~isfield(input, 'ori'); % old style
isgrad   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'coilpos')) || isgrad;             % new style
iselec   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'elecpos')) || iselec;             % new style
isnirs   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'transceiver');             
islabel  = isa(input, 'cell')    && ~isempty(input) && isa(input{1}, 'char');
haslabel = isa(input, 'struct')  && isfield(input, 'label');

if ~(isdata || isheader || isgrad || iselec || isnirs || islabel || haslabel) && isfield(input, 'hdr')
  input    = input.hdr;
  isheader = true;
end

if isdata
  % the input may be a data structure which then contains a grad/elec structure, a header or only the labels
  % preferably look at the data and not the header for the grad, because it might be re-balanced and/or planar
  if isfield(input, 'grad')
    sens   = input.grad;
    isgrad = true;
  elseif isfield(input, 'elec')
    sens   = input.elec;
    iselec = true;
  elseif issubfield(input, 'hdr.grad')
    sens   = input.hdr.grad;
    isgrad = true;
  elseif issubfield(input, 'hdr.elec')
    sens   = input.hdr.elec;
    iselec = true;
  elseif issubfield(input, 'hdr.opto')
    sens   = input.hdr.opto;
    isnirs = true;
  elseif issubfield(input, 'hdr.label')
    sens.label = input.hdr.label;
    islabel    = true;
  elseif isfield(input, 'label')
    sens.label = input.label;
    islabel    = true;
  else
    sens = [];
  end
  
elseif isheader
  if isfield(input, 'grad')
    sens   = input.grad;
    isgrad = true;
  elseif isfield(input, 'elec')
    sens   = input.elec;
    iselec = true;
  elseif isfield(input, 'opto')
    sens   = input.opto;
    isnirs = true;
  elseif isfield(input, 'label')
    sens.label = input.label;
    islabel    = true;
  end
  
elseif isgrad
  sens = input;
  
elseif iselec
  sens = input;
  
elseif isnirs
  sens = input;
  
elseif islabel
  sens.label = input;
  
elseif haslabel
  % it does not resemble anything that we had expected at this location, but it does have channel labels
  % the channel labels can be used to determine the type of sensor array
  sens.label = input.label;
  islabel    = true;
  
else
  sens = [];
end


if isfield(sens, 'type')
  % preferably the structure specifies its own type
  type = sens.type;
  
  % do not make a distinction between the neuromag data with or without space in the channel names
  if strcmp(type, 'neuromag306alt')
    type = 'neuromag306';
  elseif strcmp(type, 'neuromag122alt')
    type = 'neuromag122';
  end
  
elseif isfield(input, 'nChans') && input.nChans==1 && isfield(input, 'label') && ~isempty(regexp(input.label{1}, '^csc', 'once'))
  % this is a single channel header that was read from a Neuralynx file, might be fcdc_matbin or neuralynx_nsc
  type = 'neuralynx';
  
elseif issubfield(input, 'orig.FileHeader') &&  issubfield(input, 'orig.VarHeader')
  % this is a complete header that was read from a Plexon *.nex file using read_plexon_nex
  type = 'plexon';
  
elseif issubfield(input, 'orig.stname')
  % this is a complete header that was read from an ITAB dataset
  type = 'itab';
  
elseif issubfield(input, 'orig.sys_name')
  % this is a complete header that was read from a Yokogawa dataset
  if strcmp(input.orig.sys_name, '9ch Biomagnetometer System') || input.orig.channel_count<20
    % this is the small animal system that is installed at the UCL Ear Institute
    % see http://www.ucl.ac.uk/news/news-articles/0907/09070101
    type = 'yokogawa9';
  elseif input.orig.channel_count<160
    type = 'yokogawa64';
  elseif input.orig.channel_count<300
    type = 'yokogawa160';
  else
    % FIXME this might fail if there are many bad channels
    type = 'yokogawa440';
  end
  
elseif issubfield(input, 'orig.FILE.Ext') && strcmp(input.orig.FILE.Ext, 'edf')
  % this is a complete header that was read from an EDF or EDF+ dataset
  type = 'eeg';
  
else
  % start with unknown, then try to determine the proper type by looking at the labels
  type = 'unknown';
  
  if isgrad && isfield(sens, 'type')
    type = sens.type;
    
  elseif isgrad
    % this looks like MEG
    % revert the component balancing that was previously applied
    if isfield(sens, 'balance') && strcmp(sens.balance.current, 'comp')
      sens = undobalancing(sens);
    end
    
    % determine the type of magnetometer/gradiometer system based on the channel names alone
    % this uses a recursive call to the "islabel" section further down
    type = ft_senstype(sens.label);
    %sens_temp.type = type;
    if strcmp(type, 'unknown') %|| ~ft_senstype(sens_temp,'meg')
      % although we don't know the type, we do know that it is MEG
      type = 'meg';
    end
    
  elseif iselec
    % this looks like EEG
    
    % determine the type of eeg/acquisition system based on the channel names alone
    % this uses a recursive call to the "islabel" section further down
    type = ft_senstype(sens.label);
    %sens_temp.type = type;
    if strcmp(type, 'unknown') %|| ~ft_senstype(sens_temp,'eeg')
      % although we don't know the type, we do know that it is EEG
      type = 'eeg';
    end
    
  elseif isnirs
    % this looks like NIRS
    
    % determine the type of eeg/acquisition system based on the channel names alone
    % this uses a recursive call to the "islabel" section further down
    type = ft_senstype(sens.label);
    %sens_temp.type = type;
    if strcmp(type, 'unknown') %|| ~ft_senstype(sens_temp,'eeg')
      % although we don't know the type, we do know that it is EEG
      type = 'nirs';
    end
    
  elseif islabel
    % look only at the channel labels
    if     (mean(ismember(ft_senslabel('ant128'),         sens.label)) > 0.8)
      type = 'ant128';
    elseif (mean(ismember(ft_senslabel('ctf275'),         sens.label)) > 0.8)
      type = 'ctf275';
    elseif (mean(ismember(ft_senslabel('ctfheadloc'),     sens.label)) > 0.8)  % look at the head localization channels
      type = 'ctf275';
    elseif (mean(ismember(ft_senslabel('ctf151'),         sens.label)) > 0.8)
      type = 'ctf151';
    elseif (mean(ismember(ft_senslabel('ctf64'),          sens.label)) > 0.8)
      type = 'ctf64';
    elseif (mean(ismember(ft_senslabel('ctf275_planar'),  sens.label)) > 0.8)
      type = 'ctf275_planar';
    elseif (mean(ismember(ft_senslabel('ctf151_planar'),  sens.label)) > 0.8)
      type = 'ctf151_planar';
    elseif (mean(ismember(ft_senslabel('bti248'),         sens.label)) > 0.8) % note that it might also be a bti248grad system
      type = 'bti248';
    elseif (mean(ismember(ft_senslabel('bti148'),         sens.label)) > 0.8)
      type = 'bti148';
    elseif (mean(ismember(ft_senslabel('bti248_planar'),  sens.label)) > 0.8) % note that it might also be a bti248grad_planar system
      type = 'bti248_planar';
    elseif (mean(ismember(ft_senslabel('bti148_planar'),  sens.label)) > 0.8)
      type = 'bti148_planar';
    elseif (mean(ismember(ft_senslabel('itab28'),         sens.label)) > 0.8)
      type = 'itab28';
    elseif (mean(ismember(ft_senslabel('itab153'),        sens.label)) > 0.8)
      type = 'itab153';
    elseif (mean(ismember(ft_senslabel('itab153_planar'), sens.label)) > 0.8)
      type = 'itab153_planar';
      
      % the order is important for the different yokogawa systems, because they all share the same channel names
    elseif (mean(ismember(ft_senslabel('yokogawa440'),        sens.label)) > 0.7)
      type = 'yokogawa440';
    elseif (mean(ismember(ft_senslabel('yokogawa160'),        sens.label)) > 0.4)
      type = 'yokogawa160';
    elseif (mean(ismember(ft_senslabel('yokogawa160_planar'), sens.label)) > 0.4)
      type = 'yokogawa160_planar';
    elseif (mean(ismember(ft_senslabel('yokogawa64'),         sens.label)) > 0.4)
      type = 'yokogawa64';
    elseif (mean(ismember(ft_senslabel('yokogawa64_planar'),  sens.label)) > 0.4)
      type = 'yokogawa64_planar';
    elseif all(ismember(ft_senslabel('yokogawa9'),            sens.label))
      type = 'yokogawa9';
      
    elseif any(mean(ismember(ft_senslabel('neuromag306'),     sens.label)) > 0.4) % there are two possibilities for the channel labels: with and without a space
      type = 'neuromag306';
    elseif any(mean(ismember(ft_senslabel('neuromag122'),     sens.label)) > 0.4) % there are two possibilities for the channel labels: with and without a space
      type = 'neuromag122';
      
    elseif (mean(ismember(ft_senslabel('biosemi256'),         sens.label)) > 0.8)
      type = 'biosemi256';
    elseif (mean(ismember(ft_senslabel('biosemi128'),         sens.label)) > 0.8)
      type = 'biosemi128';
    elseif (mean(ismember(ft_senslabel('biosemi64'),          sens.label)) > 0.8)
      type = 'biosemi64';
    elseif (mean(ismember(ft_senslabel('egi256'),             sens.label)) > 0.8)
      type = 'egi256';
    elseif (mean(ismember(ft_senslabel('egi128'),             sens.label)) > 0.8)
      type = 'egi128';
    elseif (mean(ismember(ft_senslabel('egi64'),              sens.label)) > 0.8)
      type = 'egi64';
    elseif (mean(ismember(ft_senslabel('egi32'),              sens.label)) > 0.8)
      type = 'egi32';
      
      % the following check on the fraction of channels in the user's data rather than on the fraction of channels in the predefined set
    elseif (mean(ismember(sens.label, ft_senslabel('eeg1020'))) > 0.8)
      type = 'eeg1020';
    elseif (mean(ismember(sens.label, ft_senslabel('eeg1010'))) > 0.8)
      type = 'eeg1010';
    elseif (mean(ismember(sens.label, ft_senslabel('eeg1005'))) > 0.8)
      type = 'eeg1005';
      
    elseif (sum(ismember(sens.label, ft_senslabel('eeg1005'))) > 10) % Otherwise it's not even worth recognizing
      type = 'ext1020'; % this will also cover small subsets of eeg1020, eeg1010 and eeg1005
    elseif any(ismember(ft_senslabel('btiref'), sens.label))
      type = 'bti'; % it might be 148 or 248 channels
    elseif any(ismember(ft_senslabel('ctfref'), sens.label))
      type = 'ctf'; % it might be 151 or 275 channels     
    
%     elseif (mean(ismember(sens.label,    ft_senslabel('nirs'))) > 0.8)
%       type = 'nirs';
    end
  end % look at label, ori and/or pos
end % if isfield(sens, 'type')

if strcmp(type, 'unknown') && ~recursion
  % try whether only lowercase channel labels makes a difference
  if islabel
    recursion = true;
    type = ft_senstype(lower(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = lower(input.label);
    recursion = true;
    type = ft_senstype(input);
    recursion = false;
  end
end

if strcmp(type, 'unknown') && ~recursion
  % try whether only uppercase channel labels makes a difference
  if islabel
    recursion = true;
    type = ft_senstype(upper(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = upper(input.label);
    recursion = true;
    type = ft_senstype(input);
    recursion = false;
  end
end

if ~isempty(desired)
  % return a boolean flag
  switch desired
    case 'ext1020'
      type = any(strcmp(type, {'eeg1005' 'eeg1010' 'eeg1020' 'ext1020'}));
    case {'eeg' 'electrode'}
      type = any(strcmp(type, {'eeg' 'electrode' 'ant128' 'biosemi64' 'biosemi128' 'biosemi256' 'egi32' 'egi64' 'egi128' 'egi256' 'eeg1005' 'eeg1010' 'eeg1020' 'ext1020'}));
    case 'biosemi'
      type = any(strcmp(type, {'biosemi64' 'biosemi128' 'biosemi256'}));
    case 'egi'
      type = any(strcmp(type, {'egi32' 'egi64' 'egi128' 'egi256'}));
    case 'meg'
      type = any(strcmp(type, {'meg' 'magnetometer' 'ctf' 'bti' 'ctf64' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar' 'neuromag122' 'neuromag306' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar' 'bti248grad' 'bti248grad_planar' 'yokogawa9' 'yokogawa160' 'yokogawa160_planar' 'yokogawa64' 'yokogawa64_planar' 'yokogawa440' 'itab' 'itab28' 'itab153' 'itab153_planar'}));
    case 'ctf'
      type = any(strcmp(type, {'ctf' 'ctf64' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar'}));
    case 'bti'
      type = any(strcmp(type, {'bti' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar' 'bti248grad' 'bti248grad_planar'}));
    case 'neuromag'
      type = any(strcmp(type, {'neuromag122' 'neuromag306'}));
    case 'babysquid'
      type = any(strcmp(type, {'babysquid74' 'artenis123' 'magview'}));
    case 'yokogawa'
      type = any(strcmp(type, {'yokogawa160' 'yokogawa160_planar' 'yokogawa64' 'yokogawa64_planar' 'yokogawa440'}));
    case 'itab'
      type = any(strcmp(type, {'itab' 'itab28' 'itab153' 'itab153_planar'}));
    case 'meg_axial'
      % note that neuromag306 is mixed planar and axial
      type = any(strcmp(type, {'neuromag306' 'ctf64' 'ctf151' 'ctf275' 'bti148' 'bti248' 'bti248grad' 'yokogawa9' 'yokogawa64' 'yokogawa160' 'yokogawa440'}));
    case 'meg_planar'
      % note that neuromag306 is mixed planar and axial
      type = any(strcmp(type, {'neuromag122' 'neuromag306' 'ctf151_planar' 'ctf275_planar' 'bti148_planar' 'bti248_planar' 'bti248grad_planar' 'yokogawa160_planar' 'yokogawa64_planar'}));
    otherwise
      type = any(strcmp(type, desired));
  end % switch desired
end % detemine the correspondence to the desired type

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {type};
previous_argin  = current_argin;
previous_argout = current_argout;

return % ft_senstype main()
