function [hdr] = read_biosig_header(filename)

% READ_BIOSIG_HEADER reads header from EEG file using the BIOSIG
% toolbox and returns it in the FCDC framework standard format
%
% Use as
%  [hdr] = read_biosig_header(filename)
%
% The following data formats are supported: EDF, BKR, CNT, BDF, GDF,
% see for full documentation http://biosig.sourceforge.net/
%
% See also READ_BIOSIG_DATA

% Copyright (C) 2004-2012, Robert Oostenveld
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

% open the file, read the header and close it again
biosig = sopen(filename,'r');
sclose(biosig);

% the BIOSIG header is defined in full detail below
% the FCDC header should at least contain the following fields
%   hdr.Fs           sampling frequency
%   hdr.nChans       number of channels
%   hdr.nSamples     number of samples per trial
%   hdr.nSamplesPre  number of pre-trigger samples in each trial
%   hdr.nTrials      number of trials
%   hdr.label        cell-array with labels of each channel

if length(biosig.SampleRate)>1 && any(diff(biosig.SampleRate))
  error('channels with different sampling rates are not supported');
else
  hdr.Fs          = biosig.SampleRate(1);
end

if length(biosig.SPR)>1 && any(diff(biosig.SPR))
  error('channels with different number of samples are not supported');
else
  hdr.nSamples  = biosig.SPR(1);
end

hdr.nChans      = biosig.NS;
hdr.nTrials     = biosig.NRec;
hdr.nSamplesPre = 0;      % this one is not in the biosig header
hdr.label       = {};     % start with empty labels and fill them below

if hdr.nTrials>1 && hdr.nSamples==1
  % there are many trials with one sample, i.e. the data is stored in a multiplexed format instead of a block format
  % represent it in the header as a continuous format
  hdr.nSamples = hdr.nTrials;
  hdr.nTrials  = 1;
end

if isfield(biosig, 'Label')
  hdr.label = biosig.Label;
end

if length(hdr.label)~=hdr.nChans
  % make default channel labels
  hdr.label = {};
  for i=1:hdr.nChans
    hdr.label{i} = sprintf('%d', i);
  end
end

% I prefer to have them as column vector
hdr.label = hdr.label(:);

% also remember the biosig header details
hdr.orig = biosig;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the BIOSIG header always contains these elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HDR.TYPE                string          type of data format
% HDR.VERSION             string          (depends on data format)
% HDR.T0                  float[1..6]     [yyyy mm dd hh MM ss.cc] see HELP CLOCK
% HDR.NS                  integer         number of channels
% HDR.SampleRate          integer         sampling frequency in [Hz]
% HDR.NRec                integer         number of records or blocks; 1 for continous data
% HDR.SPR                 integer         samples per record
% HDR.Dur                 float           Duration (in [s]) of minimal block length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the BIOSIG header optinally contains the following elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HDR.Filter.LowPass      float           [Hz]
% HDR.Filter.HighPass     float           [Hz]
% HDR.Filter.Notch        int8            0=Off, 1=ON
% HDR.PreFilt             string          filter setting
% HDR.Label               char-array      z.B. '+C3a - C3p  '
% HDR.PhysDim             string          physical dimension e.g. 'uV'
% HDR.PhysMax             float           physical maximum
% HDR.DigMax              integer         digital maximum
% HDR.PhysMin             float           physical minimum
% HDR.DigMin              integer         digital minimum
% HDR.FLAG.TRIGGERED      int             0=no, 1=yes
% HDR.FLAG.REFERENCE      string          COM, CAR: common average reference; LOC,LAR local average ref; LAP Laplacian derivation, WGT weighted average
% HDR.Classlabel          int             0: left, 1: right, etc.
% HDR.ID.Doctor                           Identification of doctor
% HDR.ID.Hospital                         Identification of Hospital
% HDR.Patient.Name                        Name of Patient
% HDR.Patient.Age                         Age of Patient
% HDR.Patient.Sex                         Patient Gender
% HDR.Patient.Handedness                  Patient Handedness
% HDR.Patient.Medication                  Medication
% HDR.Patient.Classification              Classification of Patient
