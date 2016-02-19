function [hdr] = read_nmc_archive_k_hdr(paramfile)

% READ_NMC_ARCHIVE_K_HDR extracts 'header-data' for nmc_archive_k datasets
%
% Use as
%   hdr = read_nmc_archive_k_hdr(paramfile)
%
%
% This function specifically only reads data from one of the archived
% datasets of the Neurophysiological Mechanisms of Cognition group of
% Eric Maris, at the Donders Centre for Cognition, Radboud University,
% Nijmegen, the Netherlands. It should not be used for any other data
% format.
%
%

% Copyright (C) 2009, Roemer van der Meij
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


% Checking paramfile
if exist(paramfile,'file') ~= 2
  error('specified newparams.txt file not found');
end

% Reading samplingrate and channelnumber from paramfile
paramid = fopen(paramfile);
allfound = false;
param = [];
channelnum = []; samplingrate = [];
while allfound == false
  param = fscanf(paramid,'%s',1);
  value = fgetl(paramid);
  if strcmp(param,'samplerate')
    samplingrate = str2num(value);
  end
  if strcmp(param,'channelnum')
    channelnum = str2num(value);
  end
  if ~isempty(channelnum) && ~isempty(samplingrate)
    allfound = true;
  elseif isempty(param)
    try
      % try to find channel numbers and missing channels by dirlisting
      list = dir([paramfile(1:end-10) '***']);
      % remove any remaining dirs in the list
      list(logical([list.isdir])) = [];
      % extract names, and keep the last 3 characters
      list = {list.name};
      for iitem = 1:numel(list)
        list{iitem} = str2num(list{iitem}(end-2:end));
      end
      % remove empty items by catting to a vector
      list = [list{:}];
      
      % set channelnum and missingchan
      channelnum = list(end);
      missingchan = setdiff(1:list(end),list);
    catch
      error('error in newparams.txt, certain variables not present for subject')
    end
  end
end % while
fclose(paramid);

% Get dataformat out of paramfile and default it when not present
paramid = fopen(paramfile);
allfound = false;
param = [];
dataformat = [];
while allfound == false
  param = fscanf(paramid,'%s',1);
  value = fgetl(paramid);
  if strcmp(param,'dataformat')
    dataformat = value(3:end-1);
    if ~strcmp(dataformat, 'short') && ~strcmp(dataformat, 'int16')
      error('dataformat from newparams.txt not recognized')
    end
  end
  if ~isempty(dataformat)
    allfound = true;
  elseif isempty(param)
    dataformat = 'short'; % default-value
  end
end % while
fclose(paramid);

% Determine number of bytes per sample out of dataformat (used during data-reading)
if strcmp(dataformat, 'short') || strcmp(dataformat, 'int16')
  nBytes = 2;
else
  error('dataformat from newparams.txt not recognized')
end

% Get missing channel numbers from paramfile if not read in from above
if ~exist('missingchan','var')
  paramid = fopen(paramfile);
  allfound = false;
  param = [];
  missingchan = [];
  while allfound == false
    param = fscanf(paramid,'%s',1);
    value = fgetl(paramid);
    if strcmp(param,'missingchan')
      missingchan = strtrim(value);
    end
    if ~isempty(missingchan) || isempty(param)
      allfound = true;
    end
  end % while
  fclose(paramid);
  % Reformat missingchan for later use
  if ~isempty(missingchan)
    if ~isempty(strfind(missingchan,'-'))
      missingchan(strfind(missingchan,'-')) = ':';
    end
    missingchan = str2num(missingchan);
  end
end

% Construct channel-labels and remove missing channels
channellabels = [];
for ichan = 1:channelnum
  channellabels{ichan} = ['CH' num2str(ichan)];
end
if ~isempty(missingchan)
  channellabels(missingchan) = [];
  channelnum = length(channellabels);
end


% Determining total sample number (IMPORTANT: this assumes all channel files contain the same number of samples)
% Picking first available channel file
if length(channellabels{1} == 3)
  channelext = ['.00' channellabels{1}(3:end)];
elseif length(channellabels{1} == 5)
  channelext = ['.0' channellabels{1}(3:end)];
elseif  length(channellabels{1} == 7)
  channelext = ['.' channellabels{1}(3:end)];
end
% paramfile is either newparams.txt or .params.txt
if strcmp(paramfile(end-12:end),'newparams.txt')
  datafile = [paramfile(1:end-13) channelext];
elseif strcmp(paramfile(end-10:end),'.params.txt')
  datafile = [paramfile(1:end-11) channelext];
else
  error(['could not properly parse param-file: ' paramfile])
end
datafid = fopen(datafile,'r','l');
fseek(datafid,0,'eof');
samplesnum = ftell(datafid) / nBytes;


% Determine subject from paramfilename
slashpos = strfind(datafile, '/');
subjectname = datafile((slashpos(end-2)+1):(slashpos(end-1))-1);


% Building hdr structure
hdr = [];
hdr.Fs          = samplingrate;
hdr.nChans      = channelnum;
hdr.Sessname    = datafile((slashpos(end)+1):end-4);
hdr.label       = channellabels;
hdr.nSamples    = samplesnum;
hdr.nSamplesPre = 0;
hdr.nTrials     = 1;
hdr.dataset     = datafile;
hdr.dataformat  = dataformat;
hdr.nBytes      = nBytes;
hdr.Subject     = subjectname;
