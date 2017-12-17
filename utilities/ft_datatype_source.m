function [source] = ft_datatype_source(source, varargin)

% FT_DATATYPE_SOURCE describes the FieldTrip MATLAB structure for data that is
% represented at the source level. This is typically obtained with a beamformer of
% minimum-norm source reconstruction using FT_SOURCEANALYSIS.
%
% An example of a source structure obtained after performing DICS (a frequency
% domain beamformer scanning method) is shown here
%
%           pos: [6732x3 double]       positions at which the source activity could have been estimated
%        inside: [6732x1 logical]      boolean vector that indicates at which positions the source activity was estimated
%           dim: [xdim ydim zdim]      if the positions can be described as a 3D regular grid, this contains the
%                                       dimensionality of the 3D volume
%     cumtapcnt: [120x1 double]        information about the number of tapers per original trial
%          time: 0.100                 the latency at which the activity is estimated (in seconds)
%          freq: 30                    the frequency at which the activity is estimated (in Hz)
%           pow: [6732x120 double]     the estimated power at each source position
%     powdimord: 'pos_rpt'             defines how the numeric data has to be interpreted,
%                                       in this case 6732 dipole positions x 120 repetitions (i.e. trials)
%           cfg: [1x1 struct]          the configuration used by the function that generated this data structure
%
% Required fields:
%   - pos
%
% Optional fields:
%   - time, freq, pow, coh, eta, mom, ori, cumtapcnt, dim, transform, inside, cfg, dimord, other fields with a dimord
%
% Deprecated fields:
%   - method, outside
%
% Obsoleted fields:
%   - xgrid, ygrid, zgrid, transform, latency, frequency
%
% Historical fields:
%   - avg, cfg, cumtapcnt, df, dim, freq, frequency, inside, method,
%   outside, pos, time, trial, vol, see bug2513
%
% Revision history:
%
% (2014) The subfields in the avg and trial fields are now present in the
% main structure, e.g. source.avg.pow is now source.pow. Furthermore, the
% inside is always represented as logical vector.
%
% (2011) The source representation should always be irregular, i.e. not
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
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2013-2014, Robert Oostenveld
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

% FIXME: I am not sure whether the removal of the xgrid/ygrid/zgrid fields
% was really in 2007

% get the optional input arguments, which should be specified as key-value pairs
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest') || strcmp(version, 'upcoming')
  version = '2014';
end

if isempty(source)
  return;
end

% old data structures may use latency/frequency instead of time/freq. It is
% unclear when these were introduced and removed again, but they were never
% used by any FieldTrip function itself
if isfield(source, 'frequency'),
  source.freq = source.frequency;
  source      = rmfield(source, 'frequency');
end
if isfield(source, 'latency'),
  source.time = source.latency;
  source      = rmfield(source, 'latency');
end

switch version
  case '2014'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that it has individual source positions
    source = fixpos(source);
    
    % ensure that it is always logical
    source = fixinside(source, 'logical');
    
    % remove obsolete fields
    if isfield(source, 'method')
      source = rmfield(source, 'method');
    end
    if isfield(source, 'transform')
      source = rmfield(source, 'transform');
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

    if isfield(source, 'avg') && isstruct(source.avg) && isfield(source, 'trial') && isstruct(source.trial) && ~isempty(intersect(fieldnames(source.avg), fieldnames(source.trial)))
      % it is not possible to convert both since they have the same field names
      ft_warning('removing ''avg'', keeping ''trial''');
      source = rmfield(source, 'avg');
    end
    
    if isfield(source, 'avg') && isstruct(source.avg)
      % move the average fields to the main structure
      fn = fieldnames(source.avg);
      for i=1:length(fn)
        dat = source.avg.(fn{i});
        if isequal(size(dat), [1 size(source.pos,1)])
          source.(fn{i}) = dat';
        else
          source.(fn{i}) = dat;
        end
        clear dat
      end % j
      source = rmfield(source, 'avg');
    end
    
    if isfield(source, 'inside')
      % the inside is by definition logically indexed
      probe = find(source.inside, 1, 'first');
    else
      % just take the first source position
      probe = 1;
    end
    
    if isfield(source, 'trial') && isstruct(source.trial)
      npos = size(source.pos,1);
      
      % concatenate the fields for each trial and move them to the main structure
      fn = fieldnames(source.trial);
      
      for i=1:length(fn)
        % some fields are descriptive and hence identical over trials
        if strcmp(fn{i}, 'csdlabel')
          source.csdlabel = dat;
          continue
        end
        
        % start with the first trial
        dat    = source.trial(1).(fn{i});
        datsiz = getdimsiz(source, fn{i});
        nrpt   = datsiz(1);
        datsiz = datsiz(2:end);
        
        
        if iscell(dat)
          datsiz(1) = nrpt; % swap the size of pos with the size of rpt
          val  = cell(npos,1);
          indx = find(source.inside);
          for k=1:length(indx)
            val{indx(k)}          = nan(datsiz);
            val{indx(k)}(1,:,:,:) = dat{indx(k)};
          end
          % concatenate all data as {pos}_rpt_etc
          for j=2:nrpt
            dat = source.trial(j).(fn{i});
            for k=1:length(indx)
              val{indx(k)}(j,:,:,:) = dat{indx(k)};
            end
            
          end % for all trials
          source.(fn{i}) = val;
          
        else
          % concatenate all data as pos_rpt_etc
          val = nan([datsiz(1) nrpt datsiz(2:end)]);
          val(:,1,:,:,:) = dat(:,:,:,:);
          for j=2:length(source.trial)
            dat = source.trial(j).(fn{i});
            val(:,j,:,:,:) = dat(:,:,:,:);
          end % for all trials
          source.(fn{i}) = val;

%         else
%           siz = size(dat);
%           if prod(siz)==npos
%             siz = [npos nrpt];
%           elseif siz(1)==npos
%             siz = [npos nrpt siz(2:end)];
%           end
%           val = nan(siz);
%           % concatenate all data as pos_rpt_etc
%           val(:,1,:,:,:) = dat(:);
%           for j=2:length(source.trial)
%             dat = source.trial(j).(fn{i});
%             val(:,j,:,:,:) = dat(:);
%           end % for all trials
%           source.(fn{i}) = val;
          
        end
      end % for each field
      
      source = rmfield(source, 'trial');
      
    end % if trial
    
    % ensure that it has a dimord (or multiple for the different fields)
    source = fixdimord(source);
    
    % ensure that all data fields have the correct dimensions
    fn = getdatfield(source);
    for i=1:numel(fn)
      dimord = getdimord(source, fn{i});
      dimtok = tokenize(dimord, '_');
      dimsiz = getdimsiz(source, fn{i});
      dimsiz(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions
      if numel(dimsiz)>=3 && strcmp(dimtok{1}, 'dim1') && strcmp(dimtok{2}, 'dim2') && strcmp(dimtok{3}, 'dim3')
        % convert it from voxel-based representation to position-based representation
        try
          source.(fn{i}) = reshape(source.(fn{i}), [prod(dimsiz(1:3)) dimsiz(4:end) 1]);
        catch
          ft_warning('could not reshape %s to the expected dimensions', fn{i});
        end
      end
    end
      
    
  case '2011'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that it has individual source positions
    source = fixpos(source);
    
    % remove obsolete fields
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
    % ensure that it has individual source positions
    source = fixpos(source);
    
    % remove obsolete fields
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
    % ensure that it has individual source positions
    source = fixpos(source);
    
    % remove obsolete fields
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
    ft_error('unsupported version "%s" for source datatype', version);
end

function pos = grid2pos(xgrid, ygrid, zgrid)
[X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
pos = [X(:) Y(:) Z(:)];

function pos = dim2pos(dim, transform)
[X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
pos = [X(:) Y(:) Z(:)];
pos = ft_warp_apply(transform, pos, 'homogenous');
