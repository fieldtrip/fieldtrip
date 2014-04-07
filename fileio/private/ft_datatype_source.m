function source = ft_datatype_source(source, varargin)

% FT_DATATYPE_SOURCE describes the FieldTrip MATLAB structure for data that is
% represented at the source level. This is typically obtained with a beamformer of
% minimum-norm source reconstruction using FT_SOURCEANALYSIS.
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
%   - pos
%
% Optional fields:
%   - time, freq, pow, mom, ori, dim, cumtapcnt, cfg, dimord, other fields with a dimord
%
% Deprecated fields:
%   - method
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
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2013, Robert Oostenveld
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
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2011';
end

if isempty(source)
  return;
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
  case 'upcoming' % this is under development and expected to become the standard in 2014
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that it has individual source positions
    source = fixpos(source);
    
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
    
    if isfield(source, 'avg') && isstruct(source.avg)
      % move the average fields to the main structure
      fn = fieldnames(source.avg);
      for i=1:length(fn)
        source.(fn{i}) = source.avg.(fn{i});
      end % j
      source = rmfield(source, 'avg');
    end
    
    if isfield(source, 'trial') && isstruct(source.trial)
      npos = size(source.pos,1);
      nrpt = numel(source.trial); % note that this field is also used further doen
      
      % concatenate the fields for each trial and move them to the main structure
      fn = fieldnames(source.trial);
      
      for i=1:length(fn)
        % start with the first trial
        dum = source.trial(1).(fn{i});
        
        % some fields are descriptive and hence identical over trials
        if strcmp(fn{i}, 'csdlabel')
          source.csdlabel = dum;
          continue
        end
        
        if iscell(dum)
          val = cell(nrpt, npos);
          dum = reshape(dum, [1 npos]);
          val(1,:) = dum;
          for j=2:length(source.trial)
            dum = source.trial(j).(fn{i});
            dum = reshape(dum, [1 npos]);
            % concatenate them as {rpt_pos}
            val(j,:) = dum;
          end % for all trials
          source.(fn{i}) = val;
          
        else
          val = nan([nrpt size(dum)]);
          dum = reshape(dum, [1 size(dum)]);
          val(1,:,:,:,:) = dum;
          for j=2:length(source.trial)
            dum = source.trial(j).(fn{i});
            dum = reshape(dum, [1 size(dum)]);
            % concatenate them as rpt_pos_etc
            val(j,:,:,:,:) = dum;
          end % for all trials
          source.(fn{i}) = val;
        end
      end % for each field
      
      source = rmfield(source, 'trial');
    elseif isfield(source, 'cumtapcnt')
      % note that this field is also used further doen
      nrpt = length(source.cumtapcnt);
    else
      % note that this field is also used further doen
      nrpt = nan;
    end % if trial
    
    if isfield(source, 'cfg') && isfield(source.cfg, 'channel')
      nchan = length(source.cfg.channel);
    else
      nchan = nan;
    end
    
    % ensure that it is always logical
    source = fixinside(source, 'logical');
    
    fn = fieldnames(source);
    for i=1:length(fn)
      npos  = size(source.pos,1);
      if isfield(source, 'time')
        ntime = length(source.time);
      else
        ntime = nan;
      end
      if isfield(source, 'freq')
        nfreq = length(source.freq);
      else
        nfreq = nan;
      end
      
      dimord = [];
      
      if isnumeric(source.(fn{i}))
        val = source.(fn{i});
        switch numel(val)
          case npos
            dimord = 'pos';
          case npos*ntime
            dimord = 'pos_time';
          case nrpt*npos*ntime
            dimord = 'rpt_pos_time';
          case npos*nfreq
            dimord = 'pos_freq';
          case nrpt*npos*nfreq
            dimord = 'rpt_pos_freq';
          case npos*nfreq*ntime
            if isequal(size(val), [npos nfreq ntime])
              dimord = 'pos_freq_time';
            elseif isequal(size(val), [npos ntime nfreq])
              dimord = 'pos_time_freq';
            else
              error('cannot determine dimord for %s', fn{i});
            end
          case nrpt*npos*nfreq*ntime
            if isequal(size(val), [nrpt npos nfreq ntime])
              dimord = 'rpt_pos_freq_time';
            elseif isequal(size(val), [nrpt npos ntime nfreq])
              dimord = 'rpt_pos_time_freq';
            else
              error('cannot determine dimord for %s', fn{i});
            end
          otherwise
            dimord = [];
        end % switch
        
      elseif iscell(source.(fn{i}))
        if isfield(source, 'inside')
          % it is logically indexed
          probe = find(source.inside, 1, 'first');
        else
          % just take the first source position
          probe = 1;
        end
        
        if length(source.(fn{i}))==npos
          
          val = source.(fn{i}){probe};
          switch numel(val)
            case 1
              dimord = '{pos}';
            case ntime
              dimord = '{pos}_time';
            case 3*ntime
              dimord = '{pos}_ori_time';
            case nfreq
              dimord = '{pos}_freq';
            case 3*nfreq
              dimord = '{pos}_ori_freq';
            case 3*nrpt
              dimord = '{pos}_ori_rpt';
            case nchan*3
              dimord = '{pos}_chan_ori';
            case nfreq*ntime
              if isequal(size(val), [nfreq ntime])
                dimord = '{pos}_freq_time';
              elseif isequal(size(val), [ntime nfreq])
                dimord = '{pos}_time_freq';
              else
                error('cannot determine dimord for %s', fn{i});
              end
            case 3*nfreq*ntime
              if isequal(size(val), [3 nfreq ntime])
                dimord = '{pos}_ori_freq_time';
              elseif isequal(size(val), [3 ntime nfreq])
                dimord = '{pos}_ori_time_freq';
              else
                error('cannot determine dimord for %s', fn{i});
              end
            case 3*3
              dimord = '{pos}_ori_ori';
            otherwise
              dimord = [];
          end % switch

        else
          error('unsupported dimensions for cell-array data');
        end % if npos
        
      end % if isnumeric or iscell
      
      if ~isempty(dimord)
        source.([fn{i} 'dimord']) = dimord;
      end
    end % for each field
    
    % ensure that it has a dimord (or multiple for the different fields)
    source = fixdimord(source);
    
    
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
    error('unsupported version "%s" for source datatype', version);
end

function source = fixpos(source)
if ~isfield(source, 'pos')
  if isfield(source, 'xgrid') && isfield(source, 'ygrid') && isfield(source, 'zgrid')
    source.pos = grid2pos(source.xgrid, source.ygrid, source.zgrid);
  elseif isfield(source, 'dim') && isfield(source, 'transform')
    source.pos = dim2pos(source.dim, source.transform);
  else
    error('cannot reconstruct individual source positions');
  end
end

function pos = grid2pos(xgrid, ygrid, zgrid)
[X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
pos = [X(:) Y(:) Z(:)];

function pos = dim2pos(dim, transform)
[X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
pos = [X(:) Y(:) Z(:)];
pos = ft_warp_apply(transform, pos, 'homogenous');


