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
%
% Optional fields:
%
% Deprecated fields:
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
% $Id:$


