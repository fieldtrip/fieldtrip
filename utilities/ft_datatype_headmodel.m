function vol = ft_datatype_headmodel(vol, varargin)

% FT_DATATYPE_HEADMODEL describes the FieldTrip MATLAB structure for headmodel data
%
% The headmodel datatype represents headmodel geometry and volume conductor realated variables
% (like conductivity) typically obtained after calling FT_HEADMODEL_XXX and FT_READ_HEADSHAPE. 
% The geometry can contain one or multiple  triangulated surfaces and is normally used in the further step
% of lead field calculation in the routine FT_COMPUTE_LEADFIELD and FT_PREPARE_VOL_SENS.
%
% An example of a headshape data structure is
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
% $Id$

% get the optional input arguments, which should be specified as key-value pairs
version       = ft_getopt(varargin, 'version', 'latest');
if strcmp(version, 'latest')
  version = '2011';
end

switch version
  
  case '2011'
    if isfield(vol, 'skin_surface')
      vol.skin = vol.skin_surface;
      vol = rmfield(vol, 'skin_surface');
    end
    if strcmp(vol.type,'bem_cp')
      vol.type = 'bemcp';
    end
    if strcmp(vol.type,'openmeeg')
      vol.type = 'bem_openmeeg';
    end
    if strcmp(vol.type,'dipoli')
      vol.type = 'bem_dipoli';
    end
    if strcmp(vol.type,'asa')
      vol.type = 'bem_asa';
    end
    if strcmp(vol.type,'avo')
      vol.type = 'bem_avo';
    end
    
  case '2007'
    
  otherwise
    error('unsupported version "%s" for raw datatype', version);  
  
end
