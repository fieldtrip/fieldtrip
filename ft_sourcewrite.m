function ft_sourcewrite(cfg, source)

% FT_SOURCEWRITE exports source analysis results to gifti or nifti format,
% depending on whether the source locations are described by on a
% cortically constrained sheet (gifti) or by a regular 3D lattice (nifti)
%
% Use as
%  ft_sourcewrite(cfg, source) 
% where source is a source structure obtained after FT_SOURCEANALYSIS and 
% cfg is a configuratioun structure that should contain 
%
%  cfg.filename  = string, name of the file
%  cfg.parameter = string, functional parameter to be written to file
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this the input data will be read from a *.mat
% file on disk. This mat file should contain only a single variable, 
% corresponding with the input data structure.
%
% See also FT_SOURCEANALYSIS FT_SOURCEDESCRIPTIVES FT_VOLUMEWRITE

% Copyright (C) 2011, Jan-Mathijs Schoffelen
%
% $Id$

revision = '$Id$';

ft_defaults                 
ft_preamble init           
ft_preamble provenance        
ft_preamble trackconfig     
ft_preamble debug
ft_preamble loadvar source

source = ft_checkdata(source, 'datatype', 'source', 'feedback', 'yes');

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'parameter', 'filename'});

% get the options
cfg.parameter = ft_getopt(cfg, 'parameter');
cfg.filename  = ft_getopt(cfg, 'filename');

if isfield(source, 'dim')
  % source positions are on a regular 3D lattice, save as nifti
  if numel(cfg.filename)<=4 || ~strcmp(cfg.filename(end-3:end), '.nii');
    cfg.filename = cat(2, cfg.filename, '.nii');
  end
  
  % convert functional data into 4D
  dat = getsubfield(source, cfg.parameter);
  dat = reshape(dat, source.dim(1), source.dim(2), source.dim(3), []);
  
  if ~isfield(source, 'transform')
    source.transform = pos2transform(source.pos);
  end
  ft_write_mri(cfg.filename, dat, 'dataformat', 'nifti', 'transform', source.transform);
  
elseif isfield(source, 'tri')
  % there's a specification of a mesh, so the source positions are on a 2D
  % sheet, save as gifti
  if numel(cfg.filename)<=4 || ~strcmp(cfg.filename(end-3:end), '.gii');
    cfg.filename = cat(2, cfg.filename, '.gii');
  end
  
  if ~isfield(source, 'pnt')
    source.pnt = source.pos;
  end
  ft_write_headshape(cfg.filename, source, 'data', getsubfield(source,cfg.parameter), 'format', 'gifti');
  
else
  error('the input data does not look like a 2D sheet, nor as a 3D regular volume');
end

ft_postamble debug
ft_postamble trackconfig      
ft_postamble provenance        
