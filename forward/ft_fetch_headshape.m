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

% $Id: $

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

hdmfile   = ft_getopt(cfg,'hdmfile');
headshape = ft_getopt(cfg,'headshape');
if ischar(cfg.headshape)
  hdmfile = headshape;
end

if nargin<2
  data = [];
end

% booleans
hashdmfile    = ~isempty(hdmfile); 
hasheadshape  = ~isempty(headshape);
hasdata       = ~isempty(data);
hasdatavol    = isfield(data, 'vol');
hasdatashape  = ~isempty(data) && ~hasdatavol;

if (hashdmfile + hasheadshape + hasdata)==0
  error('No headshape found in input')
end

% check input arguments
if hasdata
  data = ft_checkdata(data);
else
  data = struct; % initialize as empty struct
end

cfg = ft_checkconfig(cfg);

if hasheadshape
  % in case headshape is given as a file string, convert it
  if ischar(headshape)
    hdmfile = headshape;
    hashdmfile = true;
    hasheadshape = false;
  elseif isa(headshape, 'config')
    headshape = struct(headshape);
  end
end

if (hashdmfile + hasheadshape + hasdata) > 1
  display = @warning;
  fprintf('Your data and/or configuration allow for multiple headmodel definitions.\n');
  keyboard
else
  display = @fprintf;
end

% get the head model definition
if hashdmfile
  display('reading headmodel from file ''%s''\n', hdmfile);
  [bnd] = ft_read_headshape(hdmfile);
elseif hasdatavol
  display('using headmodel specified in the data\n');
  vol = data.vol;
  [bnd] = vol.bnd;
elseif hasheadshape
  display('using headmodel specified in the configuration\n');
  [bnd] = getbnd(headshape);
elseif hasdatashape
  [bnd] = getbnd(data);
else
  error('no headmodel specified');
end

function bnd = getbnd(input)
bnd = [];
if isfield(input,'bnd')
  bnd = input.bnd;
elseif isfield(input,'vol') 
  bnd = input.vol.bnd;
else
  bnd = input;
end
