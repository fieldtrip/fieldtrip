function source = ft_datatype_source(source, varargin)

% FT_DATATYPE_SOURCE describes the FieldTrip MATLAB structure for source data
%
% The source datatype represents source-level data typically obtained after
% calling FT_SOURCEANALYSIS.
%
% An example of a source structure obtained after performing DICS (a frequency
% domain beamformer scanning method) is shown here
%
%           pos: [6732x3 double]       positions at which the source activity could have been estimated
%        inside: [1x3415 double]       indices to the positions at which the source activity is actually estimated
%       outside: [1x3317 double]       indices to the positions at which the source activity has not been estimated
%           dim: [xdim ydim zdim]      if the positions can be described as a 3D regular grid, this contains the
%                                       dimensionality of the 3D volume
%     cumtapcnt: [120x1 double]        information about the number of tapers per original trial
%          freq: 6                     the frequency of the oscillations at which the activity is estimated
%        method: 'singletrial'         specifies how the data is represented
%           cfg: [1x1 struct]          the configuration used by the function that generated this data structure
%           pow: [6732x120 double]     the numeric data
%     powdimord: 'pos_rpt'             defines how the numeric data has to be interpreted,
%                                       in this case 6732 dipole positions x 120 observations
%
% Required fields:
%   - pos, dimord
%
% Optional fields:
%   - time, freq, pow, mom, ori, dim, cumtapcnt, channel, other fields with a dimord
%
% Deprecated fields:
%   - method
%
% Obsoleted fields:
%   - xgrid, ygrid, zgrid, transform, latency, frequency
%
% Revision history:
%
% (2011/latest) The source representation should always be irregular, i.e. not
% a 3-D volume, contain a "pos" field and not contain a "transform".
%
% (2010) The source structure should contain a general "dimord" or specific
% dimords for each of the fields. The source reconstruction in the avg and
% trial substructures has been moved to the toplevel.
%
% (2007) The xgrid/ygrid/zgrid fields have been removed, because they are
% redundant.
%
% (2003) The initial version was defined
%
% See also FT_DATATYPE and FT_DATATYPE_xxx

% Copyright (C) 2011, Robert Oostenveld
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

% FIXME: I am not sure whether the removal of the xgrid/ygrid/zgrid fields 
% was really in 2007

% get the optional input arguments, which should be specified as key-value pairs
version = keyval('version', varargin); if isempty(version), version = 'latest'; end

if strcmp(version, 'latest')
  version = '2011';
end

% old data structures may use latency/frequency instead of time/freq. It is
% unclear when these were introduced and removed again, but they were never
% used by any fieldtrip function itself
if isfield(source, 'frequency'),
  source.freq = source.frequency;
  source      = rmfield(source, 'frequency');
end
if isfield(source, 'latency'),
  source.time = source.latency;
  source      = rmfield(source, 'latency');
end

switch version
  case '2011'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(source, 'xgrid')
      source = rmfield(source, 'xgrid');
    end
    if isfield(source, 'ygrid')
      source = rmfield(source, 'ygrid');
    end
    if isfield(source, 'zgrid')
      source = rmfield(source, 'zgrid');
    end

    if isfield(source, 'transform')
      source = rmfield(source, 'transform');
    end

    % ensure that it has a dimord (or multiple for the different fields)
    source = fixdimord(source);

  case '2010'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(source, 'xgrid')
      source = rmfield(source, 'xgrid');
    end
    if isfield(source, 'ygrid')
      source = rmfield(source, 'ygrid');
    end
    if isfield(source, 'zgrid')
      source = rmfield(source, 'zgrid');
    end

    % ensure that it has a dimord (or multiple for the different fields)
    source = fixdimord(source);

  case '2007'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(source, 'dimord')
      source = rmfield(source, 'dimord');
    end

    if isfield(source, 'xgrid')
      source = rmfield(source, 'xgrid');
    end
    if isfield(source, 'ygrid')
      source = rmfield(source, 'ygrid');
    end
    if isfield(source, 'zgrid')
      source = rmfield(source, 'zgrid');
    end

  case '2003'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(source, 'dimord')
      source = rmfield(source, 'dimord');
    end

    if ~isfield(source, 'xgrid') || ~isfield(source, 'ygrid') || ~isfield(source, 'zgrid')
      if isfield(source, 'dim')
        minx = min(source.pos(:,1));
        maxx = max(source.pos(:,1));
        miny = min(source.pos(:,2));
        maxy = max(source.pos(:,2));
        minz = min(source.pos(:,3));
        maxz = max(source.pos(:,3));
        source.xgrid = linspace(minx, maxx, source.dim(1));
        source.ygrid = linspace(miny, maxy, source.dim(2));
        source.zgrid = linspace(minz, maxz, source.dim(3));
      end
    end

  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for source datatype', version);
end

