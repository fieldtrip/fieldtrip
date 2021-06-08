function [info] = read_yorkinstruments_hdf5_meta(datafile, acq_run)

% info=read_yorkinstruments_hdf5_meta(datafile)
% Collects the required Fieldtrip header data from the data file 'filename'
%
% Copyright (C) 2020, York Neuroimaging Centre, University of York, R Aveyard
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

%hdf5 files typically contain multiple acquisitions
% including COH runs, 'default' links to the main data
if nargin < 2
  acq_run='default';
end

%Check that the datafile exists
if ~isempty(datafile),
%Check that the acquisition exists in the data file
  try
    check = h5info(datafile,strcat('/acquisitions/',num2str(acq_run)));
  catch
    error('Invalid YI HDF5 file: Missing dataset "/acquisitions/%d".',acq_run);
  end
  %Check that this acquisition is 'real' data, not COH or something else 
  info.acq_type=char(h5readatt(datafile,  strcat('/acquisitions/',num2str(acq_run))  ,'acq_type'));
  if  ~strcmp(info.acq_type,'ACQ'),
    error('"/acquisitions/%d" is not a data acquisition run.',acq_run);
  end

  info.SampleFrequency=h5readatt(datafile,  strcat('/acquisitions/',num2str(acq_run))  ,'sample_rate');
  info.Sequence=h5readatt(datafile,  strcat('/acquisitions/',num2str(acq_run))  ,'sequence');
  info.Description=char(h5readatt(datafile,  strcat('/acquisitions/',num2str(acq_run))  ,'description'));
  info.StartTime=char(h5readatt(datafile,  strcat('/acquisitions/',num2str(acq_run))  ,'start_time'));
  info.WeightsConfigured=char(h5readatt(datafile,  strcat('/acquisitions/',num2str(acq_run))  ,'weights_configured'));
  info.ChNames= h5read(datafile,[strcat('/acquisitions/',num2str(acq_run)) '/channel_list/']);
  [info.NChannels, null]=size(info.ChNames);
  Data = h5read(datafile,[strcat('/acquisitions/',num2str(acq_run)) '/data/']);
  [null, info.NSamples]=size(Data);
  try [null info.NTrials] = size(h5read(datafile,[strcat('/acquisitions/',num2str(acq_run)) '/epochs/trigger_codes']));
  catch info.NTrials=1
  end
  
  for i = 1:info.NChannels
    info.ChUnit(i)=h5readatt(datafile, char(strcat('/config/channels/',info.ChNames(i) )) ,'units');
    if info.ChUnit{i}=='?'
	    info.ChUnit{i}='unknown';
    end

    info.ChType(i)=h5readatt(datafile, char(strcat('/config/channels/',info.ChNames(i) )) ,'chan_type');
    try info.ChType(i)=h5readatt(datafile, char(strcat('/config/channels/',info.ChNames(i) )) ,'mode');
    catch
      continue
    end
  end

   
end
      





