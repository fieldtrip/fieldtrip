function ft_sourcewrite(cfg, source)

% FT_SOURCEWRITE exports source-reconstructed results to gifti or nifti format file.
% The appropriate output file depends on whether the source locations are described by
% on a cortically constrained sheet (gifti) or by a regular 3D lattice (nifti).
%
% Use as
%  ft_sourcewrite(cfg, source)
% where source is a source structure obtained from FT_SOURCEANALYSIS and
% cfg is a structure that should contain
%
%  cfg.filename  = string, filename without the extension
%  cfg.filetype  = string, can be 'nifti', 'gifti' or 'cifti' (default is automatic)
%  cfg.parameter = string, functional parameter to be written to file
%  cfg.precision = string, can be 'single', 'double', etc.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this the input data will be read from a *.mat
% file on disk. This mat file should contain only a single variable,
% corresponding with the input data structure.
%
% See also FT_SOURCEANALYSIS FT_SOURCEDESCRIPTIVES FT_VOLUMEWRITE

% Copyright (C) 2011, Jan-Mathijs Schoffelen
% Copyright (C) 2011-2014, Jan-Mathijs Schoffelen, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar source
ft_preamble provenance source
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'parameter', 'filename'});

% get the options
cfg.parameter      = ft_getopt(cfg, 'parameter');
cfg.filename       = ft_getopt(cfg, 'filename');
cfg.filetype       = ft_getopt(cfg, 'filetype');        % the default is determined further down
cfg.brainstructure = ft_getopt(cfg, 'brainstructure');  % is used for cifti
cfg.parcellation   = ft_getopt(cfg, 'parcellation');    % is used for cifti
cfg.precision      = ft_getopt(cfg, 'precision');       % is used for cifti


% check if the input data is valid for this function
if strcmp(cfg.filetype, 'cifti')

  % keep the transformation matrix
  if isfield(source, 'transform')
    transform = source.transform;
  elseif isfield(source, 'brainordinate') && isfield(source.brainordinate, 'transform')
    transform = source.brainordinate.transform;
  else
    transform = [];
  end

  if isfield(source, 'brainordinate')
    % it is a parcellated source representation, i.e. the main structure one channel for each parcel
    brainordinate = source.brainordinate;
    source = rmfield(source, 'brainordinate');

    % split them and check individually
    source        = ft_checkdata(source, 'datatype', {'timelock', 'freq', 'chan'}, 'feedback', 'yes');
    brainordinate = ft_checkdata(brainordinate, 'datatype', 'parcellation', 'parcellationstyle', 'indexed', 'hasunit', 'yes');

    % merge them again
    source = copyfields(brainordinate, source, setdiff(fieldnames(brainordinate), {'cfg'}));
  else
    source = ft_checkdata(source, 'datatype', 'source', 'hasunit', true, 'feedback', 'yes');
  end

  % keep the transformation matrix
  if ~isempty(transform)
    source.transform = transform;
  end

end % if cifti


if isempty(cfg.filetype)
  if isfield(source, 'dim')
    % source positions are on a regular 3D lattice, save as nifti
    cfg.filetype = 'nifti';
  elseif isfield(source, 'tri')
    % there is a specification of a 2D cortical sheet, save as gifti
    cfg.filetype = 'gifti';
  else
    error('the input data does not look like a 2D sheet, nor as a 3D regular volume');
  end
end

switch (cfg.filetype)
  case 'nifti'
    if numel(cfg.filename)<=4 || ~strcmp(cfg.filename(end-3:end), '.nii');
      cfg.filename = cat(2, cfg.filename, '.nii');
    end

    % convert functional data into 4D
    dat = getsubfield(source, cfg.parameter);
    dat = reshape(dat, source.dim(1), source.dim(2), source.dim(3), []);

    if ~isfield(source, 'transform')
      source.transform = pos2transform(source.pos);
    end

    % this is a bit of a kludge, but ft_write_mri can write 3D and 4D nifti
    ft_write_mri(cfg.filename, dat, 'dataformat', 'nifti', 'transform', source.transform);

  case 'gifti'
    if numel(cfg.filename)<=4 || ~strcmp(cfg.filename(end-3:end), '.gii');
      cfg.filename = cat(2, cfg.filename, '.gii');
    end

    % this is a bit of a kludge, but ft_write_headshape can write gifti including data
    ft_write_headshape(cfg.filename, source, 'data', getsubfield(source, cfg.parameter), 'format', 'gifti');

  case 'cifti'
    % brainstructure should represent the global anatomical structure, such as CortexLeft, Thalamus, etc.
    % parcellation should represent the detailled parcellation, such as BA1, BA2, BA3, etc.
    ft_write_cifti(cfg.filename, source, 'parameter', cfg.parameter, 'brainstructure', cfg.brainstructure, 'parcellation', cfg.parcellation, 'precision', cfg.precision);

  otherwise
    error('unsupported output format (%s)', cfg.filetype);
end % switch filetype

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous source
ft_postamble provenance
