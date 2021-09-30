function varargout = openpose_keypoints(filename, hdr, begsample, endsample, chanindx)

% OPENPOSE_KEYPOINTS reads keypoints from a series of OpenPose JSON files
%
% Use as
%   hdr = openpose_keypoints(filename);
%   dat = openpose_keypoints(filename, hdr, begsample, endsample, chanindx);
%   evt = openpose_keypoints(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, MOTION_C3D

% Copyright (C) 2021, Robert Oostenveld
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

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

% use the full filename including path to distinguish between similarly named files in different directories
[p, f, x] = fileparts(filename);
if isempty(p)
  % no path was specified
  fullname = which(filename);
elseif startsWith(p, ['.' filesep])
  % a relative path was specified
  fullname = fullfile(pwd, p(3:end), [f, x]);
else
  fullname = filename;
end

[p, f, x] = fileparts(fullname);

% The file names are formatted as S_Npred_02_crop_000000010636_keypoints.json
% where the long number "000000010636" corresponds to the frame number

% make a pattern for DIR with a wildcard for the frame number
tok = tokenize(f, '_');
tok{end-1} = '*';
pattern = sprintf('%s_', tok{:});     % reconstruct the filename with a wildcard and '_' between the segments
pattern = [pattern(1:end-1), x];      % remove the last '_', add the extension
dirlist  = dir(fullfile(p, pattern));
dirlist  = {dirlist.name}';

% for SSCANF we need another pattern
tok{end-1} = '%d';
pattern = sprintf('%s_', tok{:});     % reconstruct the filename with a wildcard and '_' between the segments
pattern = [pattern(1:end-1), x];      % remove the last '_', add the extension

firstframe = sscanf(dirlist{1}, pattern);
lastframe  = sscanf(dirlist{end}, pattern);
numframes  = length(dirlist);

if (lastframe-firstframe+1)~=numframes
  ft_warning('there are missing frames, the timeseries might be discontinuous');
end

% read a single file
orig = jsondecodefile(fullfile(p, dirlist{1}));

% for consistency with the file format details, "people" will be used instead of "person"
numpeople = length(orig.people);

fn = fieldnames(orig.people(1));
sel2d = contains(fn, 'keypoints_2d');
sel3d = contains(fn, 'keypoints_3d');
for field=find(sel2d)'
  dat = orig.people(1).(fn{field});
  if isempty(dat)
    numchan(field) = 0;
  else
    numchan(field) = length(dat);
  end
end
for field=find(sel3d)'
  dat = orig.people(1).(fn{field});
  if isempty(dat)
    numchan(field) = 0;
  else
    numchan(field) = length(dat);
  end
end

% construct one channel for eaach x/y/z point and for each quality metric
label = cell(sum(numpeople*numchan),1);
chantype = cell(sum(numpeople*numchan),1);
chanunit = cell(sum(numpeople*numchan),1);

indx = 1;
for people=1:numpeople
  for field=1:numel(fn)
    if sel2d(field)
      prefix = fn{field}(1:end-length('_keypoints_2d')); % the string up to _keypoints_2d
      for chan=1:numchan(field)/3
        label{indx} = sprintf('people_%d_%s_%02d_x', people, prefix, chan);
        chantype{indx} = 'position';
        chanunit{indx} = 'pixels';
        indx = indx + 1;
        label{indx} = sprintf('people_%d_%s_%02d_y', people, prefix, chan);
        chantype{indx} = 'position';
        chanunit{indx} = 'pixels';
        indx = indx + 1;
        label{indx} = sprintf('people_%d_%s_%02d_q', people, prefix, chan);
        chantype{indx} = 'quality';
        chanunit{indx} = 'unknown';
        indx = indx + 1;
      end
    elseif sel3d(field)
      prefix = fn{field}(1:end-length('_keypoints_2d')); % the string up to _keypoints_3d
      for chan=1:numchan(field)/4
        label{indx} = sprintf('people_%d_%s_%02d_x', people, prefix, chan);
        chantype{indx} = 'position';
        chanunit{indx} = 'pixels';
        indx = indx + 1;
        label{indx} = sprintf('people_%d_%s_%02d_y', people, prefix, chan);
        chantype{indx} = 'position';
        chanunit{indx} = 'pixels';
        indx = indx + 1;
        label{indx} = sprintf('people_%d_%s_%02d_z', people, prefix, chan);
        chantype{indx} = 'position';
        chanunit{indx} = 'pixels';
        indx = indx + 1;
        label{indx} = sprintf('people_%d_%s_%02d_q', people, prefix, chan);
        chantype{indx} = 'quality';
        chanunit{indx} = 'unknown';
        indx = indx + 1;      end
    end
  end % for all fields
end % for all people

% the video framerate cannot be determined from the JSON files
% FIXME there should be an option to FT_PREPROCESSING and FT_READ_HEADER to set the actual framerate
framerate = 25;
ft_warning('assuming a default framerate of 25 Hz');

if needhdr
  
  % construct the header
  hdr.Fs = framerate;
  hdr.label = label;
  hdr.chantype = chantype;
  hdr.chanunit = chanunit;
  hdr.nChans = length(label);
  hdr.nSamples = numframes;
  hdr.nSamplesPre = 0;
  hdr.nTrials = 1;
  
  % remember the details from the first JSON file
  hdr.orig = orig;
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  ft_info('reading %d frames', numframes);
  % read all files
  for frame=1:numframes
    keypoints(frame) = jsondecodefile(fullfile(p, dirlist{frame}));
  end
  
  % the rows correspond to all keypoints for all people
  % the columns correspond to frames
  dat = nan(sum(numpeople*numchan), numframes);
  
  for people=1:numpeople
    for field=1:numel(fn)
      if sel2d(field) || sel3d(field) 
        begchan = sum(numchan(1:field-1))+1;
        endchan = sum(numchan(1:field));
        chansel = begchan:endchan;
        chansel = chansel + (people-1)*sum(numchan);
        for frame=1:numframes
          if numel(keypoints(frame).people)==numpeople
            slice = keypoints(frame).people(people).(fn{field}); % this is one set of keypoints in one people for one frame
          else
            slice = [];
          end
          if ~isempty(slice)
            dat(chansel,frame) = slice;
          end
        end
      end
    end % for all frames
  end % for all people
  % return the data
  varargout = {dat};
  
elseif needevt
  ft_warning('cannot read events from %s', filename);
  
  % return the events
  varargout = {[]};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function json = jsondecodefile(filename)
fid = fopen(filename, 'rt');
str = fread(fid, [1, inf], 'char=>char');
fclose(fid);
json = jsondecode(str);

