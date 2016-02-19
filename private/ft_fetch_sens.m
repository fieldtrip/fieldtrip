function [sens] = ft_fetch_sens(cfg, data)

% FT_FETCH_SENS mimics the behaviour of FT_READ_SENS, but for a FieldTrip
% data structure or a FieldTrip configuration instead of a file on disk.
%
% Use as
%   [sens] = ft_fetch_sens(cfg, data)
% where either of the two input arguments may be empty.
%
% The positions of the sensors are specified in a gradiometer or electrode configuration or
% from a layout. The sensor configuration can be passed into this function in four ways:
%  (1) in a file whose name is passed in a configuration field, and that
%      can be imported using FT_READ_SENS,
%  (2) in a configuration field,
%  (3) in a data field, or
%  (4) in a layout file, see FT_PREPARE_LAYOUT
%
% Allowed configuration or data fields:
%   gradfile      = sensor definition file to be read for MEG data
%   elecfile      = sensor definition file to be read for EEG data
%   grad          = sensor definition from MEG data
%   elec          = sensor definition from EEG data
% 
% Allowed configuration fields:
%   layout        = reference to a layout, see FT_PREPARE_LAYOUT
%   senstype      = string, can be 'meg' or 'eeg', is used choose in combined
%                   EEG/MEG data (default='eeg')
%
% See also FT_READ_SENS, FT_PREPARE_LAYOUT, FT_FETCH_DATA

% Copyright (C) 2011, J?rn M. Horschig
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

if nargin > 1 && ~isempty(data)
  data = ft_checkdata(data);
else
  data = struct; % initialize as empty struct
end

cfg = ft_checkconfig(cfg);

% meg booleans
hasgradfile = isfield(cfg, 'gradfile'); 
hascfggrad  = isfield(cfg, 'grad');     
hasdatagrad = isfield(data, 'grad');   

% eeg booleans
haselecfile = isfield(cfg, 'elecfile'); 
hascfgelec  = isfield(cfg, 'elec');     
hasdataelec = isfield(data, 'elec');    

if hascfgelec
  haselectra = isfield(cfg.elec, 'tra') ;
end

% other
haslayout   = isfield(cfg, 'layout');
iscfgsens   = isfield(cfg, 'pnt')  || isfield(cfg, 'chanpos');
isdatasens  = isfield(data, 'pnt') || isfield(data, 'chanpos');
hassenstype = isfield(cfg, 'senstype');

if (hasgradfile || hascfggrad || hasdatagrad) && ...
    (haselecfile || hascfgelec || hasdataelec) && ~hassenstype
  error('Cannot determine whether you need gradiometer or electrode sensor definition. Specify cfg.senstype as ''MEG'' or ''EEG''');
elseif hassenstype
  switch lower(cfg.senstype)
    case 'meg'
      haselecfile = 0;
      hascfgelec  = 0;
      hasdataelec = 0;
    case 'eeg'
      hasgradfile = 0;
      hascfggrad  = 0;
      hasdatagrad = 0;
    otherwise
      error('Senstype not spicified correctly, please see the documentation of FT_FETCH_SENS');
  end;
end

if (hasgradfile + hascfggrad + hasdatagrad + ...
    haselecfile + hascfgelec + hasdataelec + ...
    haslayout + iscfgsens + isdatasens) > 1
  display = @warning;
  fprintf('Your data and configuration allow for multiple sensor definitions.\n');
else
  display = @fprintf;
end

% get the gradiometer or electrode definition
if hasgradfile
  display('reading gradiometers from file ''%s''\n', cfg.gradfile);
  sens = ft_read_sens(cfg.gradfile);
elseif hascfggrad
  display('using gradiometers specified in the configuration\n');
  sens = cfg.grad;
elseif hasdatagrad
  display('using gradiometers specified in the data\n');
  sens = data.grad;
elseif haselecfile
  display('reading electrodes from file ''%s''\n', cfg.elecfile);
  sens = ft_read_sens(cfg.elecfile);
  % only keep positions and labels in case of EEG electrodes
  dum  = sens;
  sens = [];
  if isfield(dum,'chanpos')
  sens.chanpos = dum.chanpos;
  end
  if isfield(dum,'elecpos')
  sens.elecpos = dum.elecpos;
  end
  sens.label   = dum.label;
  if isfield(dum,'unit')
    sens.unit = dum.unit;
  end
  if isfield(dum,'coordsys')
    sens.coordsys = dum.coordsys;
  end
elseif hascfgelec
  display('using electrodes specified in the configuration\n');
  sens = cfg.elec;
  % only keep positions and labels in case of EEG electrodes
  dum  = sens;
  sens = [];
  if isfield(dum,'chanpos')
  sens.chanpos = dum.chanpos;
  end
  if isfield(dum,'elecpos')
  sens.elecpos = dum.elecpos;
  end
  sens.label   = dum.label;
  if haselectra
    sens.tra = dum.tra;
  end
  if isfield(dum,'unit')
    sens.unit = dum.unit;
  end
  if isfield(dum,'coordsys')
    sens.coordsys = dum.coordsys;
  end
elseif hasdataelec
  display('using electrodes specified in the data\n');
  sens = data.elec;
  % only keep positions and labels in case of EEG electrodes
  dum  = sens;
  sens = [];
  if isfield(dum,'chanpos')
  sens.chanpos = dum.chanpos;
  end
  if isfield(dum,'elecpos')
  sens.elecpos = dum.elecpos;
  end
  sens.label   = dum.label;
  if isfield(dum,'unit')
    sens.unit = dum.unit;
  end
  if isfield(dum,'coordsys')
    sens.coordsys = dum.coordsys;
  end
elseif haslayout
  display('Using the 2-D layout to determine the sensor position\n');
  lay = ft_prepare_layout(cfg);
  sens = [];
  sens.label = lay.label;
  sens.chanpos = lay.pos;
  sens.chanpos(:,3) = 0;
elseif iscfgsens
  % could be a sensor description
  display('The configuration input might already be a sensor description.\n');
  sens = cfg;
elseif isdatasens
  % could be a sensor description
  display('The data input might already be a sensor description.\n');
  sens = data;
else
  error('no electrodes or gradiometers specified.');
end

% ensure that the sensor description is up-to-date
sens = ft_datatype_sens(sens);

