function vol = ft_datatype_headmodel(vol, varargin)

% FT_DATATYPE_HEADMODEL describes the FieldTrip MATLAB structure for a volume
% conduction model of the head that can be used for forward computations of
% the EEG potentials or the MEG fields. The volume conduction model represents
% the geometrical and the conductive properties of the head. These determine
% how the secondary (or impressed) currents flow and how these contribute to
% the model potential or field.
%
% There is a large number of forward solutions for the EEG and MEG supported
% in FieldTrip, each with its own specification of the MATLAB structure that
% describes the volume conduction model of th ehead. It would be difficult to
% list all the possibilities here. One common feature is that the volume
% conduction model should specify its type, and that preferably it should
% specify the geometrical units in which it is expressed (e.g. mm, cm or m).
%
% An example of an EEG volume conduction model with 4 councentric spheres is
% shown here:
%
% vol =
%        r: [86 88 94 100]
%        o: [0 0 0]
%     type: 'concentric'
%     unit: 'mm'
%
% For each of the methods XXX for the volume conduction model, a corresponding
% function FT_HEADMODEL_XXX exists that contains all specific details and
% references to literature that describes the implementation.
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
    if isfield(vol, 'source_surface')
      vol.source = vol.source_surface;
      vol = rmfield(vol, 'source_surface');
    end
    
    if isfield(vol, 'type')
      if strcmp(vol.type, 'bem_cp')
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
        error('This format is not supported anymore');
      end
    end
    
  otherwise
    error('converting to version "%s" is not supported', version);
end
