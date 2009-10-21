function [vol] = read_vol(filename, varargin)

% READ_VOL reads a volume conduction model from various manufacturer
% specific files. Currently supported are ASA, CTF, Neuromag, MBFYS
% and Matlab.
%
% Use as
%   vol = read_vol(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'fileformat'   string
%
% The volume conduction model is represented as a structure, and its
% contents depend on the type of model.
%
% See also TRANSFORM_VOL, PREPARE_VOL_SENS, COMPUTE_LEADFIELD

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: read_vol.m,v $
% Revision 1.4  2008/04/14 20:51:19  roboos
% fixed dependency for dipoli/ama
% added convert_units
%
% Revision 1.3  2008/03/06 09:27:54  roboos
% updated documentation
%
% Revision 1.2  2008/01/31 20:15:24  roboos
% added optional fileformat argument
%
% Revision 1.1  2008/01/28 20:10:11  roboos
% new functions based on existing fieldtrip code
%

% test whether the file exists
if ~exist(filename)
  error(sprintf('file ''%s'' does not exist', filename));
end

% get the options
fileformat = keyval('fileformat',  varargin);

% determine the filetype
if isempty(fileformat)
  fileformat = filetype(filename);
end

switch fileformat
  case 'matlab'
    matfile = filename;   % this solves a problem with the matlab compiler v3
    warning('off', 'MATLAB:load:variableNotFound');
    tmp = load(matfile, 'vol');
    warning('on', 'MATLAB:load:variableNotFound');
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
vol = convert_units(vol);
