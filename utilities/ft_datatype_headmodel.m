function vol = ft_datatype_headmodel(vol, varargin)

% FT_DATATYPE_HEADMODEL describes the FieldTrip MATLAB structure for headmodel data
%
% The headmodel datatype represents headmodel geometry and volume conductor realated variables
% (like conductivity) typically obtained after calling FT_HEADMODEL_XXX and FT_READ_HEADSHAPE. 
% The geometry can contain one or multiple  triangulated surfaces and is normally used in the further step
% of lead field calculation in the routine FT_COMPUTE_LEADFIELD and FT_PREPARE_VOL_SENS.
%
% An example of a headmodel data structure is
%
%          bnd: {struct}      the boundary (can be a struct array of nested triangulated surfaces)
%         type: 'bem_dipoli'  the type of forward solution attached to this structure
%
% Required fields:
%   - type
% 
% Optional fields:
%   - bnd, r, o, label
% 
% Deprecated fields:
%   - 
% 
% Obsoleted fields:
%
% Revision history:
% 
%
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2011, Cristiano Micheli, Robert Oostenveld
% 
% $Id: $

% get the optional input arguments, which should be specified as key-value pairs
version       = ft_getopt(varargin, 'version', 'latest');
if strcmp(version, 'latest')
  version = '2011v2';
end

if isempty(vol)
  return;
end

switch version
    
  case '2011v2' 
    if isfield(vol, 'skin_surface')
      vol.skin = vol.skin_surface;
      vol = rmfield(vol, 'skin_surface');
    end    
    if isfield(vol, 'type') && strcmp(vol.type,'bem')
      error('not able to convert the original ''bem'' volume type, try using vol.type=''dipoli'''); 
    end
    if isfield(vol, 'type') && strcmp(vol.type,'bem_cp')
      vol.type = 'bemcp';
    end
    if isfield(vol, 'type') && strcmp(vol.type,'avo')
      error('This format is not supported anymore');
    end 
    
  case {'2011v1' '2010v2' '2010v1' '2009v2' '2009v1'}
    if isfield(vol, 'skin_surface')
      vol.skin = vol.skin_surface;
      vol = rmfield(vol, 'skin_surface');
    end
    if isfield(vol, 'type') && strcmp(vol.type,'bem')
      error('not able to convert the original ''bem'' volume type, try using vol.type=''dipoli'''); 
    end
    if isfield(vol, 'type') && strcmp(vol.type,'bem_cp')
      vol.type = 'bemcp';
    end

    
  otherwise
    error('converting to version "%s" is not supported', version);  
  
end
