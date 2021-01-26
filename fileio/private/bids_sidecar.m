function sidecar = bids_sidecar(filename, suffix)

% BIDS_SIDECAR will search for corresponding BIDS sidecar files that go together with
% a specific data file. This function respects the inheritance rules and will also
% search higher up in the directory structure.
%
% Use as
%   sidecar = bids_sidecar(filename, sidecar, retval)
% where filename refers to a BIDS data file and suffix is a string that refers to the
% specific sidecar file. To read the json sidecar corresponding to the data itself,
% you can keep the suffix empty. In that case the suffix (e.g., meg or eeg) will
% be determined from the filename.
%
% This supports, but is not restricted to the following json sidecar files
%   'meg'
%   'eeg'
%   'ieeg'
%   'nirs'
%   'coordsystem'
%
% This supports, but is not restricted to the following tsv sidecar files
%   'channels'
%   'electrodes'
%   'optodes'
%   'events'
%
% In case both a tsv and a json sidecar file are present, the tsv file will be returned.
%
% See https://bids-specification.readthedocs.io/ for the specification and
% http://bids.neuroimaging.io/ for background information.
%
% See also BIDS_DATAFILE, BIDS_TSV, EVENTS_TSV, FT_READ_HEADER, FT_READ_EVENT

% Copyright (C) 2020, Robert Oostenveld
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

if nargin<2
  % automatically determine the suffix
  suffix = [];
end

if endsWith(filename, '.gz')
  % this applies to nifti files that have been gzipped
  [p, f, x] = fileparts(filename(1:end-3));
  x = [x '.gz'];
else
  [p, f, x] = fileparts(filename);
end

% get the key-value pairs that comprise the file name
entities = split(f, '_');

if isempty(suffix)
  % the suffix, i.e. the part just before the file extension, is the datatype
  suffix = entities{end};
end

% don't include the suffix in the list of entities
entities = entities(1:end-1);

if isempty(p)
  p0 = pwd;
else
  p0 = p;
end
p1 = fileparts(p0); % one directory up
p2 = fileparts(p1); % one directory up
p3 = fileparts(p2); % one directory up

if isempty(p1), p1 = fullfile(p0, '..'); end
if isempty(p2), p2 = fullfile(p1, '..'); end
if isempty(p3), p3 = fullfile(p2, '..'); end

dirlist = dir(p0);
filelist0 = {dirlist(~[dirlist.isdir]).name};
dirlist = dir(p1);
filelist1 = {dirlist(~[dirlist.isdir]).name};
dirlist = dir(p2);
filelist2 = {dirlist(~[dirlist.isdir]).name};
if numel(entities)>1 && startsWith(entities{2}, 'ses-')
  % there is a sessions level
  dirlist = dir(p3);
  filelist3 = {dirlist(~[dirlist.isdir]).name};
else
  % there is not a sessions level
  filelist3 = {};
end

% conctruct a full list of candidate files, including their full path
filelist = {};
for i=1:numel(filelist0)
  filelist{end+1} = fullfile(p0, filelist0{i});
end
for i=1:numel(filelist1)
  filelist{end+1} = fullfile(p1, filelist1{i});
end
for i=1:numel(filelist2)
  filelist{end+1} = fullfile(p2, filelist2{i});
end
for i=1:numel(filelist3)
  filelist{end+1} = fullfile(p3, filelist3{i});
end

% start with an empty return value
sidecar = [];

% if there is both a tsv and a json file for the desired sidecar, we want to return the tsv
% sort them to get the tsv files first in the list, followed by the json files
filelist = [filelist(endsWith(filelist', 'tsv')) filelist(endsWith(filelist', 'json'))];

% we are searching for a file with the datatype as suffix and that ends with json
for i=1:numel(filelist)
  [p, f, x] = fileparts(filelist{i});
  tmp = split(f, '_');
  % check the file extension, the suffix and the entities of each candidate file
  if ismember(x, {'.tsv', '.json'}) && strcmp(tmp{end}, suffix) && all(ismember(tmp(1:end-1), entities))
    ft_info('found matching BIDS sidecar ''%s''', filelist{i})
    sidecar = filelist{i};
    break % do not consider any of the other potential matches
  end
end
