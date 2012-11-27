function chanunit = ft_chanunit(input, desired)

% FT_CHANUNIT is a helper function that tries to determine the physical
% units of each channel. In case the type of channel is not detected, it
% will return 'unknown' for that channel.
%
% Use as
%   unit = ft_chanunit(hdr)
% or as
%   unit = ft_chanunit(hdr, desired)
%
% If the desired unit is not specified as second input argument, this
% function returns a Nchan*1 cell array with a string describing the
% physical units of each channel, or 'unknown' if those cannot be
% determined.
%
% If the desired unit is specified as second input argument, this function
% returns a Nchan*1 boolean vector with "true" for the channels that match
% the desired physical units and "false" for the ones that do not match.
%
% The specification of the channel units depends on the acquisition system,
% for example the neuromag306 system includes channel with the following
% units: uV, T and T/cm.
%
% See also FT_CHANTYPE

% Copyright (C) 2011, Robert Oostenveld
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


% determine the type of input, this is handled similarly as in FT_CHANTYPE
isheader =  isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'Fs');
isgrad   =  isa(input, 'struct') && isfield(input, 'pnt') && isfield(input, 'ori');
isgrad   = (isa(input, 'struct') && isfield(input, 'coilpos')) || isgrad;
isgrad   = (isa(input, 'struct') && isfield(input, 'chanpos')) || isgrad;
islabel  =  isa(input, 'cell')   && isa(input{1}, 'char');

hdr   = input;
grad  = input;
label = input;

if isheader
  numchan = hdr.nChans;
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

% start with unknown unit for all channels
chanunit = repmat({'unknown'}, size(input.label));

if ft_senstype(input, 'unknown')
  % don't bother doing all subsequent checks to determine the type of sensor array
  
elseif ft_senstype(input, 'neuromag') && isheader && issubfield(input, 'orig.chs')
  for i = 1:numchan % make a cell array of units for each channel
    switch hdr.orig.chs(i).unit
      case 201 % defined as constants by MNE, see p. 217 of MNE manual
        input.chanunit{i} = 'T/m';
      case 112
        input.chanunit{i} = 'T';
      case 107
        input.chanunit{i} = 'V';
      case 202
        input.chanunit{i} = 'Am';
      otherwise
        input.chanunit{i} = 'unknown';
    end
  end
  
elseif ft_senstype(input, 'neuromag') && isfield(input, 'chantype') && isgrad
  % look at the type of the channels
  chanunit(strcmp('eeg',              grad.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('emg',              grad.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('eog',              grad.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('ecg',              grad.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('megmag',           grad.chantype)) = {'T'};
  
  if all(sum(abs(input.tra),2)==1 | sum(abs(input.tra),2)==2)
    % it is not scaled with distance
    chanunit(strcmp('megplanar',        grad.chantype)) = {'T'};
  else
    % it is scaled with distance
    if isfield(input, 'unit')
      assumption = sprintf('T/%s', input.unit);
      chanunit(strcmp('megplanar',        input.chantype)) = {assumption};
      warning('assuming that planar channel units are %s, consistent with the geometrical units', assumption);
    else
      % the channel units remain unknown
    end
  end
  
elseif ft_senstype(input, 'neuromag') && isfield(input, 'chantype')
  % determine the units only based on the channel name and type
  chanunit(strcmp('eeg',              input.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('emg',              input.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('eog',              input.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('ecg',              input.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('megmag',           input.chantype)) = {'T'};
  if isfield(input, 'unit')
    assumption = sprintf('T/%s', input.unit);
    chanunit(strcmp('megplanar',        input.chantype)) = {assumption};
    warning('assuming that planar channel units are %s, consistent with the geometrical units', assumption);
  else
    % the channel units remain unknown
  end
  
elseif ft_senstype(input, 'ctf') && isfield(input, 'chantype')
  chanunit(strcmp('eeg',              input.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('emg',              input.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('eog',              input.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('ecg',              input.chantype)) = {'unknown'}; % FIXME
  chanunit(strcmp('meggrad',          input.chantype)) = {'T'};
  chanunit(strcmp('refmag',           input.chantype)) = {'T'};
  chanunit(strcmp('refgrad',          input.chantype)) = {'T'};
  
elseif ft_senstype(input, 'yokogawa') && isfield(input, 'chantype')
  chanunit(strcmp('meggrad',          input.chantype)) = {'unknown'}; % FIXME don't know whether it is T or T/m
  chanunit(strcmp('gegplanar',        input.chantype)) = {'unknown'}; % FIXME don't know whether it is T or T/m
  
elseif ft_senstype(input, 'bti') && isfield(input, 'chantype')
  chanunit(strcmp('meg',                 input.chantype)) = {'T'}; % this was the channel type until approx. 2 November 2012, see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1807
  chanunit(strcmp('megmag',              input.chantype)) = {'T'}; % for most 4D/BTi systems
  chanunit(strcmp('meggrad',             input.chantype)) = {'unknown'}; % FIXME don't know whether it is T or T/m
  
elseif ft_senstype(input, 'itab') && isfield(input, 'chantype')
  chanunit(strcmp('megmag',              input.chantype)) = {'T'};
  
end % if senstype

if nargin>1
  chanunit = strcmp(desired, chanunit);
end
