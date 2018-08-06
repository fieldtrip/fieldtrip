function ft_write_spike(filename, spike, varargin)

% FT_WRITE_SPIKE writes animal electrophysiology spike timestamps and/or waveforms
% to file
%
% Use as
%   ft_write_spike(filename, spike, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'dataformat'          string, see below
%   'fsample'             sampling frequency of the waveforms
%   'chanindx'            index of selected channels
%   'TimeStampPerSample'  number of timestamps per sample
%
% The supported dataformats are
%   neuralynx_nse
%   neuralynx_nts
%   plexon_nex
%   matlab
%
% See also FT_READ_SPIKE

% Copyright (C) 2007-2012, Robert Oostenveld
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

% get the options
dataformat          = ft_getopt(varargin, 'dataformat');
fsample             = ft_getopt(varargin, 'fsample');
chanindx            = ft_getopt(varargin, 'chanindx');

% FIXME rename the option TimeStampPerSample to ftimestamp, c.f. fsample
TimeStampPerSample  = keyval('TimeStampPerSample',  varargin);

% optionally select channels
if ~isempty(chanindx)
  spike.label     = spike.label(chanindx);
  spike.waveform  = spike.waveform(chanindx);
  spike.timestamp = spike.timestamp(chanindx);
end

nchans = length(spike.label);

switch dataformat
  case 'matlab'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plain MATLAB file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.mat']);
    save(filename, 'spike');

  case 'neuralynx_nts'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single channel Neuralynx NTS file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.nts']);
    if nchans>1
      ft_error('only supported for single-channel data');
    end

    label     = spike.label{1};
    timestamp = spike.timestamp{1};
    clear spike

    nts           = [];
    nts.TimeStamp = uint64(timestamp);

    % construct the file header
    nts.hdr.CheetahRev            = '4.23.0';
    % nts.hdr.NLX_Base_Class_Type   = 'SEScAcqEnt';
    nts.hdr.NLX_Base_Class_Name   = label;
    nts.hdr.RecordSize            = 8;

    % write the NSE structure to file
    write_neuralynx_nts(filename, nts);

    if 0
      % the following code snippet can be used for testing
      nts2 = read_neuralynx_nts(filename, 1, inf);
    end

  case 'neuralynx_nse'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single channel Neuralynx NSE file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.nse']);
    if nchans>1
      ft_error('only supported for single-channel data');
    end

    label     = spike.label{1};
    waveform  = spike.waveform{1};
    timestamp = spike.timestamp{1};
    unit      = spike.unit{1};
    nrecords  = length(timestamp);
    clear spike

    % cut the data into record-size pieces around the spike segments
    nrecords = length(timestamp);
    nse                   = [];
    nse.CellNumber        = unit;
    nse.TimeStamp         = uint64(timestamp);
    nse.ScNumber          = zeros(1,nrecords);
    nse.Param             = zeros(8,nrecords);
    nse.dat               = waveform;

    % construct the file header
    nse.hdr.CheetahRev            = '4.23.0';
    nse.hdr.NLX_Base_Class_Type   = 'SEScAcqEnt';
    nse.hdr.NLX_Base_Class_Name   = label;
    nse.hdr.RecordSize            = 48 + size(waveform,1)*2;  % depends on the number of samples in each waveform, normal is 112
    nse.hdr.SamplingFrequency     = fsample;

    % write the NSE structure to file
    write_neuralynx_nse(filename, nse);

    if 0
      % the following code snippet can be used for testing
      nse2 = read_neuralynx_nse(filename, 1, inf);
    end

  case 'plexon_nex'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plexon NEX file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.nex']);
    if nchans>1
      ft_error('only supported for single-channel data');
    end

    nex = [];
    nex.hdr.FileHeader.Frequency  = TimeStampPerSample*fsample;
    nex.hdr.VarHeader.Type        = 3;
    nex.hdr.VarHeader.Name        = spike.label{1};
    nex.hdr.VarHeader.WFrequency  = fsample;
    nex.var.ts  = spike.timestamp{1};
    nex.var.dat = spike.waveform{1};
    write_plexon_nex(filename, nex);

    if 0
      % the following code snippet can be used for testing
      nex2 = [];
      [nex2.var, nex2.hdr] = read_plexon_nex(filename, 'channel', 1);
    end

  otherwise
    ft_error('not implemented');
end % switch dataformat

