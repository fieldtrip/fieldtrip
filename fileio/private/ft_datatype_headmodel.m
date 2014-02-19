function vol = ft_datatype_headmodel(vol, varargin)

% FT_DATATYPE_HEADMODEL describes the FieldTrip MATLAB structure for a volume
% conduction model of the head that can be used for forward computations of
% the EEG potentials or the MEG fields. The volume conduction model represents
% the geometrical and the conductive properties of the head. These determine
% how the secondary (or impressed) currents flow and how these contribute to
% the model potential or field.
%
% A large number of forward solutions for the EEG and MEG are supported
% in FieldTrip, each with its own specification of the MATLAB structure that
% describes the volume conduction model of th ehead. It would be difficult to
% list all the possibilities here. One common feature is that the volume
% conduction model should specify its type, and that preferably it should
% specify the geometrical units in which it is expressed (e.g. mm, cm or m).
%
% An example of an EEG volume conduction model with 4 concentric spheres is:
%
% vol =
%        r: [86 88 94 100]
%        c: [0.33 1.00 0.042 0.33]
%        o: [0 0 0]
%     type: 'concentricspheres'
%     unit: 'mm'
%
% An example of an MEG volume conduction model with a single sphere fitted to
% the scalp with its center 4 cm above the line connecting the ears is:
%
% vol =
%        r: [12]
%        o: [0 0 4]
%     type: 'singlesphere'
%     unit: 'cm'
%
% For each of the methods XXX for the volume conduction model, a corresponding
% function FT_HEADMODEL_XXX exists that contains all specific details and
% references to literature that describes the implementation.
%
% Required fields:
%   - type
%
% Optional fields:
%   - unit
%
% Deprecated fields:
%   - inner_skull_surface, source_surface, skin_surface, source, skin
%
% Obsoleted fields:
%   - <none specified>
%
% Revision history:
%
% (2014/latest) All numeric values are represented in double precision.
%
% (2013) Always use the field "cond" for conductivity.
%
% (2012) Use consistent names for the volume conductor type in the structure, the
% documentation and for the actual implementation, e.g. bem_openmeeg -> openmeeg,
% fem_simbio -> simbio, concentric -> concentricspheres. Deprecated the fields
% that indicate the index of the innermost and outermost surfaces.
%
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2011-2012, Cristiano Micheli, Robert Oostenveld
%
% $Id$

% get the optional input arguments, which should be specified as key-value pairs
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2014';
end

if isempty(vol)
  return;
end

switch version
  
  case '2014'
    % first make it consistent with the 2012 version
    vol = ft_datatype_headmodel(vol, 'version', '2012');
    
    % ensure that all numbers are represented in double precision
    vol = ft_struct2double(vol);
    
  case '2013'
    % first make it consistent with the 2012 version
    vol = ft_datatype_headmodel(vol, 'version', '2012');
    
    % then rename (if neccessary the c into cond
    if isfield(vol, 'c') && ~isfield(vol, 'cond')
      vol.cond = vol.c;
      vol = rmfield(vol, 'c');
    elseif isfield(vol, 'cond') && isfield(vol, 'c') && isequal(vol.cond, vol.c)
      vol = rmfield(vol, 'c');
    elseif isfield(vol, 'cond') && isfield(vol, 'c') && ~isequal(vol.cond, vol.c)
      error('inconsistent specification of conductive properties for %s model', vol.type);
    end
    
  case '2012'
    % the following will be determined on the fly in ft_prepare_vol_sens
    if isfield(vol, 'skin_surface'),        vol = rmfield(vol, 'skin_surface');        end
    if isfield(vol, 'source_surface'),      vol = rmfield(vol, 'source_surface');      end
    if isfield(vol, 'inner_skull_surface'), vol = rmfield(vol, 'inner_skull_surface'); end
    if isfield(vol, 'skin'),                vol = rmfield(vol, 'skin');                end
    if isfield(vol, 'source'),              vol = rmfield(vol, 'source');              end
    
    % ensure a consistent naming of the volume conduction model types
    % these should match with the FT_HEADMODEL_XXX functions
    if isfield(vol, 'type')
      if strcmp(vol.type, 'concentric')
        vol.type = 'concentricspheres';
      elseif strcmp(vol.type, 'nolte')
        vol.type = 'singleshell';
      elseif strcmp(vol.type, 'multisphere')
        vol.type = 'localspheres';
      elseif strcmp(vol.type, 'bem_cp')
        vol.type = 'bemcp';
      elseif strcmp(vol.type, 'bem_dipoli')
        vol.type = 'dipoli';
      elseif strcmp(vol.type, 'bem_asa')
        vol.type = 'asa';
      elseif strcmp(vol.type, 'bem_openmeeg')
        vol.type = 'openmeeg';
      elseif strcmp(vol.type, 'fem_simbio')
        vol.type = 'simbio';
      elseif strcmp(vol.type, 'fdm_fns')
        vol.type = 'fns';
      elseif strcmp(vol.type, 'bem')
        error('not able to convert the original ''bem'' volume type, try using vol.type=''dipoli''');
      elseif strcmp(vol.type, 'avo')
        error('this format is not supported anymore');
      end
    end
    
    if isfield(vol, 'sens')
      % this applies to type=interpolate, ensure that the sensor description is up to date
      vol.sens = ft_datatype_sens(vol.sens);
    end
    
    if isfield(vol, 'type') && any(strcmp(vol.type, {'concentricspheres', 'singlesphere'}))
      if isfield(vol, 'cond') && ~isfield(vol, 'c')
        vol.c = vol.cond;
        vol = rmfield(vol, 'cond');
      elseif isfield(vol, 'cond') && isfield(vol, 'c') && isequal(vol.cond, vol.c)
        vol = rmfield(vol, 'cond');
      elseif isfield(vol, 'cond') && isfield(vol, 'c') && ~isequal(vol.cond, vol.c)
        error('inconsistent specification of conductive properties for %s model', vol.type);
      end
    end
    
    % ensure that the geometrical units are specified
    if ~isfield(vol, 'unit')
      vol = ft_convert_units(vol);
    end
    
  otherwise
    error('converting to version "%s" is not supported', version);
end
