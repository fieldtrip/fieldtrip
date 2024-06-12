function [datafile, jsonfile] = bids_datafile(filename)

% BIDS_DATAFILE will search for the data file, given one of the corresponding sidecar files
%
% Use as
%   [datafile, jsonfile] = bids-datafile(filename)
%
% See also BIDS_SIDECAR, BIDS_TSV, EVENTS_TSV

% Copyright (C) 2020-2024, Robert Oostenveld
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

if endsWith(filename, '.gz')
  % this applies to nifti files that have been gzipped
  [p, f, x] = fileparts(filename(1:end-3));
  x = [x '.gz'];
else
  [p, f, x] = fileparts(filename);
end

% get the key-value pairs that comprise the file name
entities = split(f, '_');

% the suffix is the part of the BIDS filename just before the file extension
suffix = entities{end};

% don't include the suffix in the list of entities
entities = entities(1:end-1);

% this is the list of modalities with time-series data
% FIXME in principle this function could also be made to work for other modalities
modality = {'audio', 'bold', 'eeg', 'emg', 'exg', 'eyetracker', 'ieeg', 'meg', 'motion', 'nirs', 'physio' 'stim', 'video'};

% identify the datafile that goes with this events.tsv and its corresponding json file
% this is the reverse procedure than the one implemented in BIDS_SIDECAR
for i=1:numel(modality)
  f = [sprintf('%s_', entities{:}), modality{i}];

  % we don't know the file extension, hence we have to search with a wildcard
  d = dir(fullfile(p, [f '.*']));
  if ~isempty(d)
    datafile = fullfile(p, d.name);

    % the json file can have fewer entities and/or can be at a higher level due to the inheritance principle
    % we search for the one that is the "closest by"
    foundjson = false;

    if ~foundjson
      % find the corresponding json file in the same directory as the data
      alloptions = binoall(numel(entities));
      for j=1:size(alloptions,1)
        jsonfile = fullfile(p, [sprintf('%s_', entities{alloptions(j,:)}), modality{i} '.json']);
        if exist(jsonfile, 'file')
          foundjson = true;
          break
        end
      end
    end % if ~foundjson
    if ~foundjson
      % find the corresponding json file one directory up
      alloptions = binoall(numel(entities));
      p = fileparts(p); % remove the last directory
      for j=1:size(alloptions,1)
        jsonfile = fullfile(p, [sprintf('%s_', entities{alloptions(j,:)}), modality{i} '.json']);
        if exist(jsonfile, 'file')
          foundjson = true;
          break
        end
      end
    end % if ~foundjson
    if ~foundjson
      % find the corresponding json file two directories up
      alloptions = binoall(numel(entities));
      p = fileparts(p); % remove the last directory
      for j=1:size(alloptions,1)
        jsonfile = fullfile(p, [sprintf('%s_', entities{alloptions(j,:)}), modality{i} '.json']);
        if exist(jsonfile, 'file')
          foundjson = true;
          break
        end
      end
    end % if ~foundjson
    if ~foundjson
      % find the corresponding json file three directories up
      alloptions = binoall(numel(entities));
      p = fileparts(p); % remove the last directory
      for j=1:size(alloptions,1)
        jsonfile = fullfile(p, [sprintf('%s_', entities{alloptions(j,:)}), modality{i} '.json']);
        if exist(jsonfile, 'file')
          foundjson = true;
          break
        end
      end
    end % if ~foundjson

    if ~foundjson
      error('cannot find JSON sidecar for %s', datafile);
    end
  end % ~isempty(d)
end % for each modality

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r =  binoall(n)
% this returns an 2^n by n matrix with all possible binomial options
% See also BINORND

r = zeros(2^n, n, 'logical');
for i=0:(2^n-1)
  r(i+1,:) = bitget(uint32(i), 1:n);
end
r = fliplr(r);
r = flipud(r);

