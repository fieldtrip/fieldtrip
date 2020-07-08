function out = read_nervus_data(nrvHdr, segment, range, chIdx)
% read_nervus_data  Returns data from Nicolet file.
%
%   OUT = read_nervus_data(NRVHDR, SEGMENT, RANGE, CHIDX) returns data in an n x m array of
%   doubles where n is the number of datapoints and m is the number
%   of channels.
%
%   NRVHDR is a header from the function read_nervus_header
%   SEGMENT is the segment number in the file to read from
%   RANGE is a 1x2 array with the [StartIndex EndIndex]  - default: all
%   and CHIDX is a vector of channel indeces - default: all
%
%   FILENAME is the file name of a file in the Natus/Nicolet/Nervus(TM)
%   format (originally designed by Taugagreining HF in Iceland)
%
%   Based on ieeg-portal/Nicolet-Reader
%   at https://github.com/ieeg-portal/Nicolet-Reader
%
% Copyright (C) 2016, Jan Brogger and Joost Wagenaar 
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%  FieldTrip is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  FieldTrip is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: $


if nargin == 0
  ft_error('Missing argument');
elseif nargin == 1
  segment = 1;
  range = [1 nrvHdr.targetSampleCount];
  chIdx = nrvHdr.matchingChannels(:);
elseif nargin == 2
  range = [1 nrvHdr.targetSampleCount];
  chIdx = nrvHdr.matchingChannels(:);
elseif nargin == 3
  chIdx = nrvHdr.matchingChannels(:);
end

assert(length(range) == 2, 'Range is [firstIndex lastIndex]');
assert(range(1) > 0, 'Range is must start at one or higher');
assert(length(segment) == 1, 'Segment must be single value.');

% Get cumulative sum segments.
cSumSegments = [0 cumsum([nrvHdr.Segments.duration])];

% Reopen .e file.
h = fopen_or_error(nrvHdr.filename,'r','ieee-le');

% Find sectionID for channels
lChIdx = length(chIdx);
sectionIdx = zeros(lChIdx,1);
for i = 1:lChIdx
  tmp = find(strcmp(num2str(chIdx(i)-1),{nrvHdr.StaticPackets.tag}),1);
  sectionIdx(i) = nrvHdr.StaticPackets(tmp).index;
end

% Iterate over all requested channels and populate array.
out = zeros(range(2) - range(1) + 1, lChIdx);
for i = 1 : lChIdx  
  %disp(['Reading channel ' num2str(chIdx) ' segment ' num2str(segment)]);
  % Get sampling rate for current channel
  curSF = nrvHdr.Segments(segment).samplingRate(chIdx(i));
  mult = nrvHdr.Segments(segment).scale(chIdx(i));
  
  % Find all sections
  allSectionIdx = nrvHdr.allIndexIDs == sectionIdx(i);
  allSections = find(allSectionIdx);
  
  % Find relevant sections
  sectionLengths = [nrvHdr.MainIndex(allSections).sectionL]./2;
  cSectionLengths = [0 cumsum(sectionLengths)];
  
  skipValues = cSumSegments(segment) * curSF;
  firstSectionForSegment = find(cSectionLengths > skipValues, 1) - 1 ;
  lastSectionForSegment = [];
  %lastSectionForSegment = firstSectionForSegment + ...
  %  find(cSectionLengths > curSF*nrvHdr.Segments(segment).duration,1) - 2 ;
  
  if isempty(lastSectionForSegment)
    lastSectionForSegment = length(cSectionLengths);
  end
  
  offsetSectionLengths = cSectionLengths - cSectionLengths(firstSectionForSegment);
  
  firstSection = find(offsetSectionLengths < range(1) ,1,'last');
  
  samplesInChannel = nrvHdr.Segments(segment).sampleCount;
  if range(2) > samplesInChannel
    endRange = samplesInChannel;
  else
    endRange = range(2);
  end
  
  lastSection = find(offsetSectionLengths >= endRange,1)-1;
  
  if isempty(lastSection)
    lastSection = length(offsetSectionLengths)-1;
  end
  
  if lastSection > lastSectionForSegment
    ft_error('Index out of range for current section: %i > %i, on channel: %i', ...
      range(2), cSectionLengths(lastSectionForSegment+1), chIdx(i));
  end
  
  useSections = allSections(firstSection: lastSection);
  useSectionL = sectionLengths(firstSection: lastSection);
  
  % First Partial Segment
  curIdx = 1;
  curSec = nrvHdr.MainIndex(useSections(1));
  fseek(h, curSec.offset,'bof');
  
  firstOffset = range(1) - offsetSectionLengths(firstSection);
  lastOffset = min([range(2) useSectionL(1)]);
  lsec = lastOffset-firstOffset + 1;
  
  fseek(h, (firstOffset-1) * 2,'cof');
  out(1 : lsec,i) = fread(h, lsec, 'int16') * mult;
  curIdx = curIdx +  lsec;
  
  if length(useSections) > 1
    % Full Segments
    for j = 2: (length(useSections)-1)
      curSec = nrvHdr.MainIndex(useSections(j));
      fseek(h, curSec.offset,'bof');
      
      out(curIdx : (curIdx + useSectionL(j) - 1),i) = ...
        fread(h, useSectionL(j), 'int16') * mult;
      curIdx = curIdx +  useSectionL(j);
    end
    
    % Final Partial Segment
    curSec = nrvHdr.MainIndex(useSections(end));
    fseek(h, curSec.offset,'bof');
    %Sometimes there is less data available in the file.
    %So we calculate the wanted lastReadLength,
    %but assign in the data array the actualReadLength
    lastReadLength = length(out)-curIdx+1;
    %disp(['LastReadLength: ' num2str(lastReadLength)]);
    lastData = fread(h, lastReadLength, 'int16');
    actualReadLength = length(lastData);
    %if (actualReadLength <lastReadLength)
    %  warning('Some data not found');
    %disp(['ActualLastReadLength: ' num2str(lastReadLength)]);    
    out(curIdx : (curIdx+actualReadLength-1),i) = lastData * mult;
  end
  
end

% Close the .e file.
fclose(h);

end


