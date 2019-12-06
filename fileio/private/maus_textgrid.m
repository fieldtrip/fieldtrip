function [varargout] = maus_textgrid(filename, hdr, begsample, endsample, chanindx)

% MAUS_TEXTGRID reads speech segments from a file that has been processed with MAUS
% see https://clarin.phonetik.uni-muenchen.de/BASWebServices
%
% Use as
%   hdr = maus_textgrid(filename);
%   dat = maus_textgrid(filename, hdr, begsample, endsample, chanindx);
%   evt = maus_textgrid(filename, hdr);
%
% You should pass the *.TextGrid file as the filename, There should be a
% corresponding wav file with the same filename in the same directory.
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV

% Copyright (C) 2019 Robert Oostenveld
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

[p, f, x] = fileparts(filename);
wavfile = fullfile(p, [f '.wav']);

if needhdr
  %% return the audio file header
  varargout{1} = ft_read_header(wavfile);
  
elseif needdat
  %% return the audio file data
  varargout{1} = ft_read_data(wavfile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);
  
elseif needevt
  %% return the segmentations in the textgrid file as events
  
  % File type = "ooTextFile"
  % Object class = "TextGrid"
  %
  % xmin = 0
  % xmax = 47.982585
  % tiers? <exists>
  % size = 3
  % item []:
  %     item [1]:
  %         class = "IntervalTier"
  %         name = "ORT-MAU"
  %         xmin = 0
  %         xmax = 47.982585
  %         intervals: size = 120
  %         intervals [1]:
  %             xmin = 0.000000
  %             xmax = 3.570000
  %             text = ""
  %         intervals [2]:
  %             xmin = 3.570000
  %             xmax = 3.730000
  %             text = "Could"
  %         intervals [3]:
  %             xmin = 3.730000
  %             xmax = 3.820000
  %             text = "he"
  %         intervals [4]:
  %             xmin = 3.820000
  %             xmax = 4.030000
  %             text = "not"
  %         intervals [5]:
  %             xmin = 4.030000
  %             xmax = 4.330000
  %             text = "now"
  %         intervals [6]:
  %             xmin = 4.330000
  %             xmax = 4.690000
  %             text = "turn"
  %         intervals [7]:
  %             xmin = 4.690000
  %             xmax = 5.160000
  %             text = "back"
  %         intervals [8]:
  %             xmin = 5.160000
  %             xmax = 5.620000
  %             text = ""
  %         intervals [9]:
  %             xmin = 5.620000
  %             xmax = 6.120000
  %             text = "Acknowledge"
  %         intervals [10]:
  %             xmin = 6.120000
  %             xmax = 6.330000
  %             text = "his"
  %         intervals [11]:
  %             xmin = 6.330000
  %             xmax = 6.620000
  %             text = "error"
  %  ...
  
  
  fid = fopen_or_error(filename, 'rt');
  
  lines = {};
  while ~feof(fid)
    lines = cat(1, lines, fgetl(fid));
  end
  fclose(fid);
  
  % remove empty lines
  lines = lines(~cellfun(@isempty, lines));
  % remove whitespace
  lines = cellfun(@strtrim, lines, 'UniformOutput', false);
  
  Nitem = str2double(getitem(lines, 'size'));
  
  for i=1:Nitem
    [dum, pos] = getitem(lines, sprintf('item [%d]:', i));
    lines = lines(pos:end);
    item(i).class = getitem(lines, 'class');
    item(i).name = getitem(lines, 'name');
    item(i).xmin = str2double(getitem(lines, 'xmin'));
    item(i).xmax = str2double(getitem(lines, 'xmax'));
    item(i).size = str2double(getitem(lines, 'intervals: size'));
    for j=1:item(i).size
      [dum, pos] = getitem(lines, sprintf('intervals [%d]:', j));
      section = lines(pos:pos+3);
      item(i).intervals(j).xmin = str2double(getitem(section, 'xmin'));
      item(i).intervals(j).xmax = str2double(getitem(section, 'xmax'));
      item(i).intervals(j).text = getitem(section, 'text');
    end
  end
  
  k = 1;
  event = [];
  for i=1:numel(item)
    for j=1:numel(item(i).intervals)
      event(k).type     = [item(i).class ' ' item(i).name];
      event(k).value    = item(i).intervals(j).text;
      event(k).sample   = round(item(i).intervals(j).xmin * hdr.Fs) + 1; % the first sample is 1
      event(k).duration = round((item(i).intervals(j).xmax - item(i).intervals(j).xmin) * hdr.Fs);
      event(k).offset   = 0;
      k = k+1;
    end
  end
  
  varargout{1} = event;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val, pos] = getitem(lines, key)
pos = find(startsWith(lines, key), 1, 'first');
pieces = split(lines{pos}, '=');
if numel(pieces)==1
  val = [];
elseif numel(pieces)==2
  val = strtrim(pieces{2});
  val(val=='"') = [];
else
  ft_error('cannot parse string');
end

