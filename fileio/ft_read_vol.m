function [vol] = ft_read_vol(filename, varargin)

% FT_READ_VOL reads a volume conduction model from various manufacturer
% specific files. Currently supported are ASA, CTF, Neuromag, MBFYS
% and Matlab.
%
% Use as
%   vol = ft_read_vol(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'fileformat'   string
%
% The volume conduction model is represented as a structure, and its
% contents depend on the type of model.
%
% See also FT_TRANSFORM_VOL, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2008-2010 Robert Oostenveld
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

% test whether the file exists
if ~exist(filename)
  error(sprintf('file ''%s'' does not exist', filename));
end

% get the options
fileformat = keyval('fileformat',  varargin);

% determine the filetype
if isempty(fileformat)
  fileformat = ft_filetype(filename);
end

switch fileformat
  case 'matlab'
    matfile = filename;   % this solves a problem with the matlab compiler v3
    ws = warning('off', 'MATLAB:load:variableNotFound');
    tmp = load(matfile, 'vol');
    warning(ws);
    vol = getfield(tmp, 'vol');

  case 'ctf_hdm'
    vol = read_ctf_hdm(filename);

  case 'asa_vol'
    vol = read_asa_vol(filename);
    vol.type = 'asa';

  case 'mbfys_ama'
    ama = loadama(filename);
    vol = ama2vol(ama);

  case 'neuromag_fif'
    % do not read the volume into Matlab, but use external Neuromag toolbox
    vol.type     = 'neuromag';
    vol.filename = filename;
    vol.chansel  = [];  % this is defined later based on the channels present in the data
    % initialize the Neuromag toolbox, this requires a gradfile and hdmfile
    fprintf('using Neuromag volume conductor from %s\n', filename);
    fprintf('using Neuromag gradiometer definition from %s\n', cfg.gradfile);
    megmodel('head', cfg.gradfile, filename);
    % read the triangulated boundary from the neuromag BEM model
    [vol.bnd.pnt, vol.bnd.tri, vol.bnd.nrm] = loadtri(vol.filename);
    vol.bnd.pnt = vol.bnd.pnt*100;  % convert to cm

  otherwise
    error('unknown fileformat for volume conductor model');
end

% this will add the units to the volume conductor model
vol = ft_convert_units(vol);
