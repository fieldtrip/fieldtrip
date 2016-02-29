function [nsout] = read_neuroshare(filename, varargin)

% READ_NEUROSHARE reads header information or data from any file format
% supported by Neuroshare. The file can contain event timestamps, spike
% timestamps and waveforms, and continuous (analog) variable data.
%
% Use as:
%   hdr = read_neuroshare(filename, ...)
%   dat = read_neuroshare(filename, ...)
%
% Optional input arguments should be specified in key-value pairs and may include:
%   'dataformat'    = string
%   'readevent'     = 'yes' or 'no' (default)
%   'readspike'     = 'yes' or 'no' (default)
%   'readanalog'    = 'yes' or 'no' (default)
%   'chanindx'      = list with channel indices to read
%   'begsample      = first sample to read
%   'endsample      = last sample to read
%
% NEUROSHARE: http://www.neuroshare.org is a site created to support the
% collaborative development of open library and data file format
% specifications for neurophysiology and distribute open-source data
% handling software tools for neuroscientists.
%
% Note that this is a test version, WINDOWS only

% Copyright (C) 2009-2011, Saskia Haegens
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

% check the availability of the required neuroshare toolbox
ft_hastoolbox('neuroshare', 1);

% get the optional input arguments
dataformat    = ft_getopt(varargin, 'dataformat');
begsample     = ft_getopt(varargin, 'begsample');
endsample     = ft_getopt(varargin, 'endsample');
chanindx      = ft_getopt(varargin, 'chanindx');
readevent     = ft_getopt(varargin, 'readevent', 'no');
readspike     = ft_getopt(varargin, 'readspike', 'no');
readanalog    = ft_getopt(varargin, 'readanalog', 'no');

% determine the filetype
if isempty(dataformat)
  if filetype_check_extension(filename, '.nev') % to prevent confusion with neuralynx nev files
    dataformat = 'nev';
  else
    dataformat = ft_filetype(filename);
  end
end

% determine the required neuroshare library
switch dataformat
  case 'nev'
    lib = 'nsNEVLibrary.dll';
  case 'plexon_plx'
    lib = 'nsPlxLibrary.dll';
  case 'xxx'
    lib = 'xxx.dll';
end



% NEUROSHARE LIBRARY %
% set the library (note: library has to be on the path)
[feedback] = ns_SetLibrary(which(lib));
if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end

% feedback about datafile and library
[feedback, libinfo] = ns_GetLibraryInfo;
if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
disp(['Loading file ' filename ' using ' libinfo.Description ' from ' libinfo.Creator]);

% open dataset
[feedback fileID] = ns_OpenFile(filename);
if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end



% HEADER %
% retrieve dataset information
[feedback hdr.fileinfo] = ns_GetFileInfo(fileID);
if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end

% retrieve entity information
[feedback hdr.entityinfo] = ns_GetEntityInfo(fileID, 1:hdr.fileinfo.EntityCount);
if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end

% hdr.entityinfo.EntityType specifies the type of entity data recorded on that
% channel. It can be one of the following:
%   0: Unknown entity
%   1: Event entity
%   2: Analog entity
%   3: Segment entity
%   4: Neural event entity
enttyp = {'unknown'; 'event'; 'analog'; 'segment'; 'neural'};
for i=1:5
  list.(enttyp{i}) = find([hdr.entityinfo.EntityType] == i-1); % gives channel numbers for each entity type
end

% give warning if unkown entities are found
if ~isempty(list.unknown)
  warning(['There are ' num2str(length(list.unknown)) ' unknown entities found, these will be ignored.']);
end

% retrieve event information
if ~isempty(list.event)
  [feedback hdr.eventinfo] = ns_GetEventInfo(fileID, list.event);
  if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
end

% retrieve analog information
if ~isempty(list.analog)
  [feedback hdr.analoginfo] = ns_GetAnalogInfo(fileID, list.analog);
  if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
end

% retrieve segment information
if ~isempty(list.segment)
  [feedback hdr.seginfo] = ns_GetSegmentInfo(fileID, list.segment);
  if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
  [feedback hdr.segsourceinfo] = ns_GetSegmentSourceInfo(fileID, list.segment, 1);
  if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
end

% retrieve neural information
if ~isempty(list.neural)
  [feedback hdr.neuralinfo] = ns_GetNeuralInfo(fileID, list.neural);
  if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
end


% required to get actual analog chan numbers
if ~isempty(list.analog)
  [feedback analog.contcount] = ns_GetAnalogData(fileID, list.analog, 1, max([hdr.entityinfo(list.analog).ItemCount]));
end


% EVENT %
% retrieve events
if strcmp(readevent, 'yes') && ~isempty(list.event)
  [feedback event.timestamp event.data event.datasize] = ns_GetEventData(fileID, list.event, 1:max([hdr.entityinfo(list.event).ItemCount]));
  if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
   
  % skip empty ones ???
  event.timestamp(event.datasize==0)=[];
  event.data(event.datasize==0)=[];
  event.datasize(event.datasize==0)=[];

  for c=1:length(list.event)
    if hdr.entityinfo(list.event(c)).ItemCount~=0
      for i=1:length(event.timestamp)
        [feedback, event.sample(i,c)] = ns_GetIndexByTime(fileID, list.event(c), event.timestamp(i,c), 0);
        if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
      end
    end
  end
elseif strcmp(readevent, 'yes') && isempty(list.event)
  warning('no events were found in the data')
end



% DATA %
% retrieve analog data
if strcmp(readanalog, 'yes') && ~isempty(list.analog)
  % set defaults
  if isempty(chanindx)
    chanindx = list.analog(analog.contcount~=0);
  else % rebuild chanindx such that doesnt contain empty channels
    if length(chanindx)==length(list.analog(analog.contcount~=0)) && chanindx(1)==1
      chanindx = list.analog(analog.contcount~=0); % only read nonempty channels
    end
    if length(chanindx)>length(list.analog(analog.contcount~=0))
      chanindx = list.analog(chanindx & analog.contcount~=0); % only read nonempty channels
    end
  end
  if isempty(begsample); 
    begsample = 1;
  end
  if isempty(endsample);
    itemcount = max([hdr.entityinfo(list.analog).ItemCount]);
  else
    itemcount = endsample - begsample + 1;
  end
  [feedback analog.contcount analog.data] = ns_GetAnalogData(fileID, chanindx, begsample, itemcount);
  if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
elseif strcmp(readanalog, 'yes') && isempty(list.analog)
  warning('no analog events were found in the data')
end



% SPIKE %
% retrieve segments       [ (sorted) spike waveforms ]
if strcmp(readspike, 'yes') && ~isempty(list.segment)
  % collect data: all chans (=list.segment) and all waveforms (=hdr.entityinfo.ItemCount)
  [feedback segment.timestamp segment.data segment.samplecount segment.unitID] = ns_GetSegmentData(fileID, list.segment, 1:max([hdr.entityinfo(list.segment).ItemCount]));
  if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
elseif strcmp(readspike, 'yes') && isempty(list.segment)
  warning('no spike waveforms were found in the data')
end

% retrieve neural events  [ spike timestamps ]
if strcmp(readspike, 'yes') && ~isempty(list.neural)
  neural.data = nan(length(list.neural), max([hdr.entityinfo(list.neural).ItemCount])); % pre-allocate
  for chan=1:length(list.neural) % get timestamps
    [feedback neural.data(chan,1:hdr.entityinfo(list.neural(chan)).ItemCount)] = ns_GetNeuralData(fileID, list.neural(chan), 1, hdr.entityinfo(list.neural(chan)).ItemCount);
    if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end
  end
elseif strcmp(readspike, 'yes') && isempty(list.neural)
  warning('no spike timestamps were found in the data')
end



% close dataset
[feedback] = ns_CloseFile(fileID);
if feedback~=0, [feedback err] = ns_GetLastErrorMsg; disp(err), end


% collect the output
nsout = [];
nsout.hdr  = hdr;
nsout.list = list;
try nsout.event  = event;    end
try nsout.analog = analog;   end
try nsout.spikew = segment;  end
try nsout.spiket = neural;   end




