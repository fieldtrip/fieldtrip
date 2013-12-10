function [dat] = read_neuralynx_sdma(dataset, begsample, endsample, chanindx);

% READ_NEURALYNX_SDMA read specified channels and samples from a Neuralynx splitted DMA dataset
%
% Use as
%    [hdr] = read_neuralynx_sdma(dataset)
%    [dat] = read_neuralynx_sdma(dataset, begsample, endsample, chanindx)
%
% The splitted DMA dataset is not a formal Neuralynx format, but at
% the FCDC we use it in conjunction with SPIKEDOWNSAMPLE. The dataset
% directory contains files, one for each channel, each containing a
% 8-byte header followed by the binary values for all samples. Commonly
% the binary values are represented as int32, but it is possible to use
% int16 or other numeric representations. The 8-byte header specifies the
% numeric representation and the bitshift that should be applied (in case
% of integer representations).
%
% This function returns the electrophysiological data in AD units
% and not in uV. You should look up the details of the headstage and
% the Neuralynx amplifier and scale the values accordingly.

% Copyright (C) 2006-2008, Robert Oostenveld
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

needhdr = (nargin==1);
needdat = (nargin>=2);

if ~isdir(dataset)
  error('dataset should be a directory');
end

% determine the content of the dataset directory
dirlist = dir(dataset);

% determine the correspondence between files and channels
label         = cell(length(dirlist),1);
filesel       = zeros(length(dirlist),1);
headerfile    = [];
[pd, nd, xd]  = fileparts(dataset);
for i=1:length(dirlist)
  [pf, nf, x1] = fileparts(dirlist(i).name);
  [pf, nf, x2] = fileparts(nf);
  % compare the name of the channel with the name of the dataset and look at the file extention
  % the file should be called "dataset/dataset.chanlabel" or "dataset/dataset.chanlabel.bin"
  if strcmp(nf, nd) && ~isempty(x1)
    if ~isempty(x2)
      label{i}   = x2(2:end);
      filesel(i) = 1;
    elseif strcmp(x1, '.txt')
      % this is not a proper channel, but an ascii file with header info
      label{i}   = [];
      filesel(i) = 0;
      headerfile = fullfile(dataset, dirlist(i).name);
    else      
      label{i}   = x1(2:end);
      filesel(i) = 1;
    end
  else
    label{i}   = [];
    filesel(i) = 0;
  end
end
label = label(filesel==1);
bytes = cell2mat({dirlist(filesel==1).bytes}); % number of bytes in each file
filelist = {dirlist(filesel==1).name};         % filename excluding the path
clear dirlist filesel pd nd xd

for i=1:length(filelist)
  % include the full path in the filename
  filelist{i} = fullfile(dataset, filelist{i});
end

% there is a preferred ordering, which corresponds with the channel order in the non-splitted DMA file
preferred = {
  'stx'
  'pid'
  'siz'
  'tsh'
  'tsl'
  'cpu'
  'ttl'
  'x01'
  'x02'
  'x03'
  'x04'
  'x05'
  'x06'
  'x07'
  'x08'
  'x09'
  'x10'
  'csc001'
  'csc002'
  'csc003'
  'csc004'
  'csc005'
  'csc006'
  'csc007'
  'csc008'
  'csc009'
  'csc010'
  'csc011'
  'csc012'
  'csc013'
  'csc014'
  'csc015'
  'csc016'
  'csc017'
  'csc018'
  'csc019'
  'csc020'
  'csc021'
  'csc022'
  'csc023'
  'csc024'
  'csc025'
  'csc026'
  'csc027'
  'csc028'
  'csc029'
  'csc030'
  'csc031'
  'csc032'
  'csc033'
  'csc034'
  'csc035'
  'csc036'
  'csc037'
  'csc038'
  'csc039'
  'csc040'
  'csc041'
  'csc042'
  'csc043'
  'csc044'
  'csc045'
  'csc046'
  'csc047'
  'csc048'
  'csc049'
  'csc050'
  'csc051'
  'csc052'
  'csc053'
  'csc054'
  'csc055'
  'csc056'
  'csc057'
  'csc058'
  'csc059'
  'csc060'
  'csc061'
  'csc062'
  'csc063'
  'csc064'
  'csc065'
  'csc066'
  'csc067'
  'csc068'
  'csc069'
  'csc070'
  'csc071'
  'csc072'
  'csc073'
  'csc074'
  'csc075'
  'csc076'
  'csc077'
  'csc078'
  'csc079'
  'csc080'
  'csc081'
  'csc082'
  'csc083'
  'csc084'
  'csc085'
  'csc086'
  'csc087'
  'csc088'
  'csc089'
  'csc090'
  'csc091'
  'csc092'
  'csc093'
  'csc094'
  'csc095'
  'csc096'
  'csc097'
  'csc098'
  'csc099'
  'csc100'
  'csc101'
  'csc102'
  'csc103'
  'csc104'
  'csc105'
  'csc106'
  'csc107'
  'csc108'
  'csc109'
  'csc110'
  'csc111'
  'csc112'
  'csc113'
  'csc114'
  'csc115'
  'csc116'
  'csc117'
  'csc118'
  'csc119'
  'csc120'
  'csc121'
  'csc122'
  'csc123'
  'csc124'
  'csc125'
  'csc126'
  'csc127'
  'csc128'
  'csc129'
  'csc130'
  'csc131'
  'csc132'
  'csc133'
  'csc134'
  'csc135'
  'csc136'
  'csc137'
  'csc138'
  'csc139'
  'csc140'
  'csc141'
  'csc142'
  'csc143'
  'csc144'
  'csc145'
  'csc146'
  'csc147'
  'csc148'
  'csc149'
  'csc150'
  'csc151'
  'csc152'
  'csc153'
  'csc154'
  'csc155'
  'csc156'
  'csc157'
  'csc158'
  'csc159'
  'csc160'
  'csc161'
  'csc162'
  'csc163'
  'csc164'
  'csc165'
  'csc166'
  'csc167'
  'csc168'
  'csc169'
  'csc170'
  'csc171'
  'csc172'
  'csc173'
  'csc174'
  'csc175'
  'csc176'
  'csc177'
  'csc178'
  'csc179'
  'csc180'
  'csc181'
  'csc182'
  'csc183'
  'csc184'
  'csc185'
  'csc186'
  'csc187'
  'csc188'
  'csc189'
  'csc190'
  'csc191'
  'csc192'
  'csc193'
  'csc194'
  'csc195'
  'csc196'
  'csc197'
  'csc198'
  'csc199'
  'csc200'
  'csc201'
  'csc202'
  'csc203'
  'csc204'
  'csc205'
  'csc206'
  'csc207'
  'csc208'
  'csc209'
  'csc210'
  'csc211'
  'csc212'
  'csc213'
  'csc214'
  'csc215'
  'csc216'
  'csc217'
  'csc218'
  'csc219'
  'csc220'
  'csc221'
  'csc222'
  'csc223'
  'csc224'
  'csc225'
  'csc226'
  'csc227'
  'csc228'
  'csc229'
  'csc230'
  'csc231'
  'csc232'
  'csc233'
  'csc234'
  'csc235'
  'csc236'
  'csc237'
  'csc238'
  'csc239'
  'csc240'
  'csc241'
  'csc242'
  'csc243'
  'csc244'
  'csc245'
  'csc246'
  'csc247'
  'csc248'
  'csc249'
  'csc250'
  'csc251'
  'csc252'
  'csc253'
  'csc254'
  'csc255'
  'csc256'
  'crc'
  };

[sel1, sel2] = match_str(preferred, label);
if length(sel2)==length(label)
  % all reorder the labels/files
  label = label(sel2);
  bytes = bytes(sel2);
  filelist = filelist(sel2);
else
  % not all files in this SDMA dataset could be accounted for in the
  % preference list, hence do not attempt to apply the preferred ordering
end

if needhdr
  % construct a general header for the complete dataset
  if ~isempty(headerfile)
    orig             = neuralynx_getheader(headerfile);
    hdr.Fs           = orig.SamplingFrequency;
    hdr.nSamples     = [];        % see below
    hdr.nSamplesPre  = 0;         % number of pre-trigger samples in each trial
    hdr.nTrials      = 1;         % number of trials
    hdr.label        = label;
    hdr.nChans       = length(label);
    hdr.orig.dataset = orig;       % keep the header details
  else
    % some parts of the header have to be hardcoded, since the splitted dataset does not contain all header information
    hdr.Fs           = 32556;     % sampling frequency
    hdr.nSamples     = [];        % see below
    hdr.nSamplesPre  = 0;         % number of pre-trigger samples in each trial
    hdr.nTrials      = 1;         % number of trials
    hdr.label        = label;
    hdr.nChans       = length(label);
  end
  
  % read the header of each individual file, i.e. for each variable
  clear orig
  for i=1:length(filelist)
    orig(i) = read_neuralynx_bin(filelist{i});
  end

  % determine the number of samples for each channel
  nsamples = cell2mat({orig.nSamples});
  if any(nsamples~=nsamples(1))
    error('different number of samples over channels are not supported');
  else
    hdr.nSamples = nsamples(1);
  end

  % determine the sampling frequency for each channel
  fsample = cell2mat({orig.Fs});
  if any(diff(fsample))
    error('different sampling rates over channels are not supported');
  elseif any(fsample~=hdr.Fs)
    error('inconsistent sampling rates');
  end
  
  % determine the first and last timestamp, by reading them from the timestamp channels
  tslfile = filelist{find(strcmp('tsl', label))};
  tshfile = filelist{find(strcmp('tsh', label))};
  beg_tsl = read_neuralynx_bin(tslfile, 1, 1);
  beg_tsh = read_neuralynx_bin(tshfile, 1, 1);
  end_tsl = read_neuralynx_bin(tslfile, hdr.nSamples, hdr.nSamples);
  end_tsh = read_neuralynx_bin(tshfile, hdr.nSamples, hdr.nSamples);
  hdr.FirstTimeStamp = timestamp_neuralynx(beg_tsl, beg_tsh);
  hdr.LastTimeStamp  = timestamp_neuralynx(end_tsl, end_tsh);
  hdr.TimeStampPerSample = double(hdr.LastTimeStamp-hdr.FirstTimeStamp)./(hdr.nSamples-1);  % this should be double, since it can be fractional

  % also remember the original header details
  hdr.orig.chan = orig;
  
  % only return the header information
  dat = hdr;

else
  % allocate memory for all data
  dat = zeros(length(chanindx), endsample-begsample+1);

  % read all channels, one small chunk at at time, and write it to separate files
  for i=1:length(chanindx)
    j = chanindx(i);
    dat(i,:) = double(read_neuralynx_bin(filelist{j}, begsample, endsample));
  end
end
