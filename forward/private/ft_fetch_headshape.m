function [bnd] = ft_fetch_headshape(cfg, data)
% FT_FETCH_HEADSHAPE reads headshape definitions from a FieldTrip
% data structure or a FieldTrip configuration instead of a file on disk.
%
% Use as
%   [bnd] = ft_fetch_headshape(cfg, data)
%
% Either of the two input arguments may be empty.
%
% The headshape definitions are specified in the configuration or
% data. The data can be passed into this function in three ways:
%  (1) in a file whose name is passed in a configuration field.
%  (2) in a configuration field,
%  (3) in a data field.
% 
% You should specify the headshape with
%   cfg.hdmfile       = string, file containing the headshape
% or alternatively
%   cfg.headshape     = structure with headshape
%   data.vol          = data with volume conduction model
%   data.bnd          = data with headshape
%
%
% See also FT_PREPARE_VOL_SENS, FT_READ_VOL, FT_FETCH_VOL

% Copyright (C) 2012, Cristiano Micheli
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

% check input arguments
if nargin > 1 && ~isempty(data)
  data = ft_checkdata(data);
else
  data = struct; % initialize as empty struct
end

cfg = ft_checkconfig(cfg);

% in case headshape is given as a file string, convert it
if ischar(cfg.headshape)
  cfg.hdmfile = cfg.headshape;
  cfg=rmfield(cfg,'headshape');
end

% booleans
hashdmfile       = isfield(cfg, 'hdmfile'); 
hascfgheadshape  = isfield(cfg, 'headshape');
hasdatavol       = isfield(data, 'vol');
hasdataheadshape = isfield(data, 'bnd');

if hashdmfile
  if isempty(cfg.hdmfile)
    cfg = rmfield(cfg,'hdmfile');
    hashdmfile = false;
  end
end

if hascfgheadshape
  if isempty(cfg.headshape)
    cfg = rmfield(cfg,'headshape');
    hascfgheadshape = false;
  else
    if isa(cfg.headshape, 'config')
      % convert the nested config-object back into a normal structure
      cfg.headshape = struct(cfg.headshape);
    end
  end
end

if (hashdmfile + hascfgheadshape + hasdatavol + hasdataheadshape) > 1
  display = @warning;
  fprintf('Your data and configuration allow for multiple headmodel definitions.\n');
  keyboard
else
  display = @fprintf;
end

% get the head model definition
if isfield(cfg, 'hdmfile')
  display('reading headmodel from file ''%s''\n', cfg.hdmfile);
  [bnd] = ft_read_headshape(cfg.hdmfile);
elseif isfield(cfg, 'headshape')
  display('using headmodel specified in the configuration\n');
  [bnd] = cfg.headshape;
elseif isfield(data, 'vol')
  display('using headmodel specified in the data\n');
  vol = data.vol;
  [bnd] = vol.bnd;
elseif isfield(data, 'bnd')
  [bnd] = data.bnd;
else
  error('no headmodel specified');
end
