function [headmodel] = ft_datatype_headmodel(headmodel, varargin)

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
% headmodel =
%        r: [86 88 94 100]
%        c: [0.33 1.00 0.042 0.33]
%        o: [0 0 0]
%     type: 'concentricspheres'
%     unit: 'mm'
%
% An example of an MEG volume conduction model with a single sphere fitted to
% the scalp with its center 4 cm above the line connecting the ears is:
%
% headmodel =
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
% (2015/latest) Use the field name "pos" instead of "pnt" for vertex positions.
%
% (2014) All numeric values are represented in double precision.
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

% get the optional input arguments, which should be specified as key-value pairs
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2015';
end

if isempty(headmodel)
  return;
end

if iscell(headmodel)
  % this might represent combined EEG, ECoG and/or MEG
  for i=1:numel(headmodel)
    % call recursively
    headmodel{i} = ft_datatype_headmodel(headmodel{i}, varargin{:});
  end
  return
end

switch version

  case '2015'
    % first make it consistent with the 2014 version
    headmodel = ft_datatype_headmodel(headmodel, 'version', '2013');

    % rename pnt into pos
    headmodel = fixpos(headmodel);

  case '2014'
    % first make it consistent with the 2013 version
    headmodel = ft_datatype_headmodel(headmodel, 'version', '2013');

    % ensure that all numbers are represented in double precision
    headmodel = ft_struct2double(headmodel);

  case '2013'
    % first make it consistent with the 2012 version
    headmodel = ft_datatype_headmodel(headmodel, 'version', '2012');

    % then rename (if neccessary the c into cond
    if isfield(headmodel, 'c') && ~isfield(headmodel, 'cond')
      headmodel.cond = headmodel.c;
      headmodel = rmfield(headmodel, 'c');
    elseif isfield(headmodel, 'cond') && isfield(headmodel, 'c') && isequal(headmodel.cond, headmodel.c)
      headmodel = rmfield(headmodel, 'c');
    elseif isfield(headmodel, 'cond') && isfield(headmodel, 'c') && ~isequal(headmodel.cond, headmodel.c)
      error('inconsistent specification of conductive properties for %s model', headmodel.type);
    end

  case '2012'
    % the following will be determined on the fly in ft_prepare_vol_sens
    if isfield(headmodel, 'skin_surface'),        headmodel = rmfield(headmodel, 'skin_surface');        end
    if isfield(headmodel, 'source_surface'),      headmodel = rmfield(headmodel, 'source_surface');      end
    if isfield(headmodel, 'inner_skull_surface'), headmodel = rmfield(headmodel, 'inner_skull_surface'); end
    if isfield(headmodel, 'skin'),                headmodel = rmfield(headmodel, 'skin');                end
    if isfield(headmodel, 'source'),              headmodel = rmfield(headmodel, 'source');              end

    % ensure a consistent naming of the volume conduction model types
    % these should match with the FT_HEADMODEL_XXX functions
    if isfield(headmodel, 'type')
      if strcmp(headmodel.type, 'concentric')
        headmodel.type = 'concentricspheres';
      elseif strcmp(headmodel.type, 'nolte')
        headmodel.type = 'singleshell';
      elseif strcmp(headmodel.type, 'multisphere')
        headmodel.type = 'localspheres';
      elseif strcmp(headmodel.type, 'bem_cp')
        headmodel.type = 'bemcp';
      elseif strcmp(headmodel.type, 'bem_dipoli')
        headmodel.type = 'dipoli';
      elseif strcmp(headmodel.type, 'bem_asa')
        headmodel.type = 'asa';
      elseif strcmp(headmodel.type, 'bem_openmeeg')
        headmodel.type = 'openmeeg';
      elseif strcmp(headmodel.type, 'fem_simbio')
        headmodel.type = 'simbio';
      elseif strcmp(headmodel.type, 'fdm_fns')
        headmodel.type = 'fns';
      elseif strcmp(headmodel.type, 'bem')
        error('not able to convert the original ''bem'' volume type, try using headmodel.type=''dipoli''');
      elseif strcmp(headmodel.type, 'avo')
        error('this format is not supported anymore');
      end
    end

    if isfield(headmodel, 'sens')
      % this applies to type=interpolate, ensure that the sensor description is up to date
      headmodel.sens = ft_datatype_sens(headmodel.sens);
    end

    if isfield(headmodel, 'type') && any(strcmp(headmodel.type, {'concentricspheres', 'singlesphere'}))
      if isfield(headmodel, 'cond') && ~isfield(headmodel, 'c')
        headmodel.c = headmodel.cond;
        headmodel = rmfield(headmodel, 'cond');
      elseif isfield(headmodel, 'cond') && isfield(headmodel, 'c') && isequal(headmodel.cond, headmodel.c)
        headmodel = rmfield(headmodel, 'cond');
      elseif isfield(headmodel, 'cond') && isfield(headmodel, 'c') && ~isequal(headmodel.cond, headmodel.c)
        error('inconsistent specification of conductive properties for %s model', headmodel.type);
      end
    end

    % ensure that the geometrical units are specified
    if ~isfield(headmodel, 'unit')
      headmodel = ft_convert_units(headmodel);
    end

  otherwise
    error('converting to version "%s" is not supported', version);
end
