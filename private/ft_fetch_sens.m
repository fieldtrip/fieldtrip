function [sens] = ft_fetch_sens(cfg, data)

% FT_FETCH_SENS mimics the behavior of FT_READ_SENS, but for a FieldTrip
% data structure or a FieldTrip configuration instead of a file on disk.
%
% Use as
%   [sens] = ft_fetch_sens(cfg)
% or as
%   [sens] = ft_fetch_sens(cfg, data)
%
% The sensor configuration can be passed into this function in four ways:
%  (1) in a configuration field
%  (2) in a file whose name is passed in a configuration field, see FT_READ_SENS
%  (3) in a layout file, see FT_PREPARE_LAYOUT
%  (4) in a data field
%
% The following fields are used from the configuration:
%   cfg.elec     = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad     = structure with gradiometer definition or filename, see FT_READ_SENS
%   cfg.opto     = structure with optode definition or filename, see FT_READ_SENS
%   cfg.layout   = structure with layout definition or filename, see FT_PREPARE_LAYOUT
%   cfg.senstype = string, can be 'meg', 'eeg', or 'nirs', this is used to choose in combined data (default = 'eeg')
%
% When the sensors are not specified in the configuration, this function will
% fetch the grad, elec or opto field from the data.
%
% See also FT_READ_SENS, FT_DATATYPE_SENS, FT_FETCH_DATA, FT_PREPARE_LAYOUT

% Copyright (C) 2011-2016, Jorn M. Horschig
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
%    but WITHOUT ANY WARRANTY; without evft_neighbourploten the implied warranty of
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed', {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed', {'optofile', 'opto'});

% set the defaults
cfg.senstype = ft_getopt(cfg, 'senstype');

% meg booleans
hasgradfile = isfield(cfg, 'grad') && ischar(cfg.grad);
hascfggrad  = isfield(cfg, 'grad') && isstruct(cfg.grad);
hasdatagrad = isfield(data, 'grad');

% eeg booleans
haselecfile = isfield(cfg, 'elec') && ischar(cfg.elec);
hascfgelec  = isfield(cfg, 'elec') && isstruct(cfg.elec);
hasdataelec = isfield(data, 'elec');

% nirs booleans
hasoptofile = isfield(cfg, 'opto') && ischar(cfg.opto);
hascfgopto  = isfield(cfg, 'opto') && isstruct(cfg.opto);
hasdataopto = isfield(data, 'opto');

% other
haslayout   = isfield(cfg, 'layout');
iscfgsens   = isfield(cfg, 'pnt')  || isfield(cfg, 'chanpos');
isdatasens  = isfield(data, 'pnt') || isfield(data, 'chanpos');

if isempty(cfg.senstype) && ((hasgradfile || hascfggrad || hasdatagrad) + (haselecfile || hascfgelec || hasdataelec) + (hasoptofile || hascfgopto || hasdataopto))>1
  ft_error('Cannot determine which sensors you want to work on. Specify cfg.senstype as ''meg'', ''eeg'' or ''nirs''');

elseif ~isempty(cfg.senstype)
  if iscell(cfg.senstype)
    % this represents combined EEG and MEG sensors, where each modality has its own sensor definition
    % use recursion to fetch all sensor descriptions
    sens = cell(size(cfg.senstype));
    for i=1:numel(cfg.senstype)
      tmpcfg = cfg;
      tmpcfg.senstype = cfg.senstype{i};
      sens{i} = ft_fetch_sens(tmpcfg, data);
    end
    return

  else
    switch lower(cfg.senstype)
      case 'meg'
        haselecfile = false;
        hascfgelec  = false;
        hasdataelec = false;
        hasoptofile = false;
        hascfgopto  = false;
        hasdataopto = false;
      case 'eeg'
        hasgradfile = false;
        hascfggrad  = false;
        hasdatagrad = false;
        hasoptofile = false;
        hascfgopto  = false;
        hasdataopto = false;
      case 'nirs'
        haselecfile = false;
        hascfgelec  = false;
        hasdataelec = false;
        hasgradfile = false;
        hascfggrad  = false;
        hasdatagrad = false;
      otherwise
        ft_error('unsupported specification of cfg.senstype as "%s"', cfg.senstype);
    end
  end
end

if (hasgradfile + hascfggrad + hasdatagrad + ...
    haselecfile + hascfgelec + hasdataelec + ...
    hasoptofile + hascfgopto + hasdataopto + ...
    haslayout + iscfgsens + isdatasens) > 1
  fprintf('Your data and configuration allow for multiple sensor definitions.\n');
  display = @warning;
else
  display = @fprintf;
end

if hasgradfile
  display('reading gradiometers from file ''%s''\n', cfg.grad);
  sens = ft_read_sens(cfg.grad, 'senstype', 'meg');
elseif hascfggrad
  display('using gradiometers specified in the configuration\n');
  sens = cfg.grad;
elseif hasdatagrad
  display('using gradiometers specified in the data\n');
  sens = data.grad;

elseif haselecfile
  display('reading electrodes from file ''%s''\n', cfg.elec);
  sens = ft_read_sens(cfg.elec, 'senstype', 'eeg');
  % only keep positions and labels in case of EEG electrodes
  sens = keepfields(sens, {'elecpos', 'chanpos', 'unit', 'coordsys', 'label','tra'});
elseif hascfgelec
  display('using electrodes specified in the configuration\n');
  sens = cfg.elec;
  % only keep positions and labels in case of EEG electrodes
  sens = keepfields(sens, {'elecpos', 'chanpos', 'unit', 'coordsys', 'label','tra'});
elseif hasdataelec
  display('using electrodes specified in the data\n');
  sens = data.elec;
  % only keep positions and labels in case of EEG electrodes
  sens = keepfields(sens, {'elecpos', 'chanpos', 'unit', 'coordsys', 'label','tra'});

elseif hasoptofile
  display('reading optodes from file ''%s''\n', cfg.opto);
  sens = ft_read_sens(cfg.opto, 'senstype', 'nirs');
  % only keep known fields in case of NIRS optodes
  sens = keepfields(sens, {'optopos', 'optotype', 'chanpos', 'unit', 'coordsys', 'label', 'transceiver', 'wavelength', 'transmits', 'laserstrength'});
elseif hascfgopto
  display('using optodes specified in the configuration\n');
  sens = cfg.opto;
  % only keep known fields in case of NIRS optodes
  sens = keepfields(sens, {'optopos', 'optotype', 'chanpos', 'unit', 'coordsys', 'label', 'transceiver', 'wavelength', 'transmits', 'laserstrength'});
elseif hasdataopto
  display('using optodes specified in the data\n');
  sens = data.opto;
  % only keep known fields in case of NIRS optodes
  sens = keepfields(sens, {'optopos', 'optotype', 'chanpos', 'unit', 'coordsys', 'label', 'transceiver', 'wavelength', 'transmits', 'laserstrength'});

elseif haslayout
  display('Using the 2-D layout to determine the sensor position\n');
  lay = ft_prepare_layout(cfg);

  % remove the COMNT and SCALE labels
  sel = ~ismember(lay.label, {'COMNT' 'SCALE'});

  sens = [];
  sens.label = lay.label(sel);
  sens.chanpos = lay.pos(sel,:);
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
  ft_error('no electrodes, gradiometers or optodes specified.');
end

% ensure that the sensor description is up-to-date
if (hasgradfile + hascfggrad + hasdatagrad + ...
    haselecfile + hascfgelec + hasdataelec + ...
    hasoptofile + hascfgopto + hasdataopto)
  % this should only be called if the sensor definition is a complete one, and not constructed from a layout
  % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3143#c9
  sens = ft_datatype_sens(sens);
end
