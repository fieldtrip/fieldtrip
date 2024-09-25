function fieldtrip2ctf(filename, data, varargin)

% FIELDTRIP2CTF saves a FieldTrip data structure to a CTF dataset.
%
% The file to which the data is exported depends on the input data structure that you
% provide. The "raw" and "timelock" structures can be exported to a CTF dataset. The
% "montage" structure can be exported to a CTF "Virtual Channels" file.
%
% Use as
%   fieldtrip2ctf(filename, data, ...)
% where filename is a string and data is a FieldTrip raw, timelock or montage
% structure.
%
% Additional options should be specified in key-value pairs and can be
%   'ds' = struct, original dataset information as obtained with readCTFds
%
% See also FT_DATATYPE, FT_APPLY_MONTAGE, FT_VOLUMEWRITE, FT_SOURCEWRITE, FT_WRITE_DATA

% Copyright (C) 2015-2024, Robert Oostenveld
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

% this ensures that the path is correct and that the ft_defaults global variable is available
ft_defaults

ds = ft_getopt(varargin, 'ds', []);

type = ft_datatype(data);

if isempty(ds) && isfield(data, 'hdr') && isfield(data.hdr, 'orig') && isfield(data.hdr.orig, 'res4')
  % take it from the data structure
  ds = data.hdr.orig;
end

switch type
  case {'raw', 'timelock'}
    % these are saved to a *.ds dataset

    % this only works if the channel names are exactly the same
    if length(data.label)~=size(ds.res4.chanNames,1)
      error('channel count is not consistent with the original dataset')
    end
    for i=1:numel(data.label)
      assert(startsWith(deblank(ds.res4.chanNames(i,:)), data.label{i}), 'channel name is not consistent with the original dataset');
    end

    if strcmp(type, 'timelock')
      ds.res4.no_trials = 1;
      if isfield(data, 'dof')
        n = unique(data.dof(:));
        if length(n)==1
          ds.res4.no_trials_avgd = n;
        else
          ft_warning('setting "no_trials_avgd" to 1')
        end
      else
        ds.res4.no_trials_avgd = 1;
      end
      % convert timelock data to raw
      data = ft_checkdata(data, 'datatype', 'raw');
    else
      ds.res4.no_trials = length(data.trial);
      ds.res4.no_trials_avgd = 0;
    end

    % concatenate the trials along the 3rd dimension
    dat = cat(3, data.trial{:});
    % permute to get samples*channels*trials
    dat = permute(dat, [2 1 3]);

    assert(~isempty(ds));

    [p, f, x] = fileparts(filename);

    % some of the original dataset fields (probably) do not apply any more
    ds.baseName     = f;
    ds.path         = p;
    ds.mrk          = [];
    ds.TrialClass   = [];
    ds.badSegments  = [];
    ds.BadChannels  = [];
    ds.processing   = [];
    ds.newds        = [];

    % some of the original dataset fields need to be updated
    ds.res4.no_samples   = size(dat,1);
    ds.res4.no_channels  = length(data.label);
    ds.res4.sample_rate  = data.fsample;

    if exist(filename, 'dir')
      ft_warning('overwriting existing dataset %s', filename);
      rmdir(filename, 's');
    end

    writeCTFds(filename, ds, dat, 'T');

  case 'montage'
    % virtual channels are weighted linear combinations of real channels
    % collected by the CTF MEG System

    % rename it for convenience
    montage = data; clear data;

    fid = fopen(filename, 'wt');
    assert(fid>0, 'could not open file for writing');
    fprintf(fid, '// Virtual channel configuration\n');
    for i=1:length(montage.labelnew)
      sel = find(montage.tra(i,:) ~= 0);
      fprintf(fid, '\n');
      fprintf(fid, 'VirtualChannel\n');
      fprintf(fid, '{\n');
      fprintf(fid, 'Name:\t%s\n', montage.labelnew{i});
      fprintf(fid, 'Unit:\n');
      for j=1:numel(sel)
        fprintf(fid, 'Ref:\t%s,%f\n', montage.labelold{sel(j)}, montage.tra(i,sel(j)));
      end
      fprintf(fid, '}\n');
    end

  otherwise
    ft_error('unsuported data structure "%s" for exporting to CTF', type);
end
