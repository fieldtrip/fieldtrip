function write_fcdc_spike(filename, spike, varargin);

% WRITE_FCDC_SPIKE writes spike timestamps and/or waveforms to a file
%
% Use as
%   write_fcdc_spike(filename, spike, ...)
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
% See also READ_FCDC_SPIKE

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: write_fcdc_spike.m,v $
% Revision 1.5  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.4  2007/10/08 13:01:47  roboos
% select channels based on chanindx
% updated documentation
%
% Revision 1.3  2007/03/21 17:25:54  roboos
% renamed neuralynx_nte to nts, since that is the correct extension
%
% Revision 1.2  2007/03/20 17:12:02  roboos
% added documentation
% implemented plexon_nex for spike waveforms
% moved scaling and conversion to int16 to low level function (for neuralynx)
% removed old pieces of code that was commented out
%

fieldtripdefs

% get the options
dataformat          = keyval('dataformat',          varargin);
fsample             = keyval('fsample',             varargin);
chanindx            = keyval('chanindx',            varargin);

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
    % plain matlab file
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
      error('only supported for single-channel data');
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
      error('only supported for single-channel data');
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
      error('only supported for single-channel data');
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
    error('not implemented');
end % switch dataformat

