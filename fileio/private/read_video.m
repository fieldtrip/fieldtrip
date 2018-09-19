function out = read_video(filename, varargin)

% READ_VIDEO
%
% Use as
%   hdr = read_video(filename)
% or
%   dat = read_video(filename, hdr, begsample, endsample)
%
% See also READ_VIDEOMEG_VID, LOAD_VIDEO123

% Copyright (C) 2015 Robert Oostenveld
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

needhdr = numel(varargin)==0;
needdat = numel(varargin)>0;

video = VideoReader(filename);

switch video.VideoFormat
  case 'RGB24'
    % there are three values for each pixel in each frame
    pixval = 3;
  otherwise
    ft_error('unsupported video format %s', video.VideoFormat);
end % switch

if needhdr
  % convert it into a FieldTrip header
  hdr         = [];
  hdr.Fs      = video.FrameRate;
  
  hdr.nChans  = video.Width * video.Height * pixval;
  
  % create artificial channels for each pixel
  hdr.label    = cellstr(num2str((1:hdr.nChans)'));
  hdr.chantype = repmat({'video'}  , size(hdr.label));
  hdr.chanunit = repmat({'unknown'}, size(hdr.label));
  
  if isfield(video, 'NumberOfFrames')
    hdr.nSamples  = video.NumberOfFrames;
  else
    hdr.nSamples  = floor(video.Duration*video.FrameRate);
  end
  hdr.nTrials     = 1;
  hdr.nSamplesPre = 0;
  
  % store it as a structure
  hdr.orig = struct(video);
  
  % return the header as output
  out = hdr;
  
elseif needdat
  hdr       = varargin{1};
  begsample = varargin{2};
  endsample = varargin{3};
  nsamples  = endsample - begsample + 1;
  
  if nargin>3
    chanindx = varargin{4};
    chanindx = chanindx(:); % ensure it to be a column
  else
    chanindx = [];
  end
  
  % read the first frame to determine the type (probably uint8)
  f = readFrame(video);
  
  % give some information, but not too often
  ft_info timeout 60
  ft_info once
  ft_info('the original size of individual video frames is [%d %d %d]', size(f));
  
  if ~isempty(chanindx)
    f = f(chanindx);
  else
    f = f(:); % ensure it to be a column
  end
  
  % allocate enough memory for all frames
  dat = repmat(f, 1, nsamples);
  
  % we cannot skip frames when reading
  for frame=2:endsample
    f = readFrame(video);
    if ~isempty(chanindx)
      f = f(chanindx);
    else
      f = f(:); % ensure it to be a column
    end
    dat(:,frame) = f;
  end
  dat = dat(:, begsample:endsample);
  
  % return the data as output
  out = dat;
  
end

delete(video)
