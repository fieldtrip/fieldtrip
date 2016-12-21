function [out] = read_tmri_poly5(filename, varargin)

% READ_TMSI_POLY5
%
% Use as
%   hdr = read_tmri_poly5(filename)
%   dat = read_tmsi_poly5(filename, hdr, begblock, endblock)
%
% This implementation is as closely as possible based on the original "tms_read",
% which contains the comments
%
%   Changed on 08-10-2014 by TL, TMSi
%   - Now supports loading a file from different directory than the script file
%   - Feedback on the validity of arguments and whether a file could be found or not.
%   - Dialogue is opened when no argument was given.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

needhdr = (nargin==1);
needdat = (nargin>1);

if needdat
  signal = varargin{1};
  begblock = varargin{2};
  endblock = varargin{3};
end

fid = fopen(filename, 'rb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if needhdr
  pos = 31;
  fseek(fid,pos,'bof');
  version = fread(fid,1,'int16');
  if version == 203
    frewind(fid);
    signal.header.FID                   = fread(fid,31,'uchar');
    signal.header.VersionNumber         = fread(fid, 1,'int16');
  else % version 204
    frewind(fid);
    signal.header.FID                   = fread(fid,32,'uchar');
    signal.header.VersionNumber         = fread(fid, 1,'int16');
  end
  signal.header.MeasurementName       = fread(fid,81,'uchar');
  signal.header.FS                    = fread(fid, 1,'int16');
  signal.header.StorageRate           = fread(fid, 1,'int16');
  signal.header.StorageType           = fread(fid, 1,'uchar');
  signal.header.NumberOfSignals       = fread(fid, 1,'int16');
  signal.header.NumberOfSamplePeriods = fread(fid, 1,'int32');
  signal.header.EMPTYBYTES            = fread(fid, 4,'uchar');
  signal.header.StartMeasurement      = fread(fid,14,'uchar');
  signal.header.NumberSampleBlocks    = fread(fid, 1,'int32');
  signal.header.SamplePeriodsPerBlock = fread(fid, 1,'uint16');
  signal.header.SizeSignalDataBlock   = fread(fid, 1,'uint16');
  signal.header.DeltaCompressionFlag  = fread(fid, 1,'int16');
  signal.header.TrailingZeros         = fread(fid,64,'uchar');
  
  %conversion to char of text values
  signal.header.FID               = char(signal.header.FID);
  signal.header.MeasurementName   = char(signal.header.MeasurementName(2:signal.header.MeasurementName(1)+1));
  %    signal.header.StartMeasurement  = typecast(int8(signal.header.StartMeasurement),'int16'); % commented by RS 1feb11
  
  signal.fs = signal.header.FS;
  
  %check for right fileversion
  % if signal.header.VersionNumber ~= 203
  %    disp('Wrong file version! Imput file must be a Poly5/TMS version 2.03 file!');
  %    return;
  % end;
  
  % Signal description
  for g=1:signal.header.NumberOfSignals,
    signal.description(g).SignalName        = fread(fid,41,'uchar');
    signal.description(g).Reserved          = fread(fid, 4,'uchar');
    signal.description(g).UnitName          = fread(fid,11,'uchar');
    signal.description(g).UnitLow           = fread(fid, 1,'float32');
    signal.description(g).UnitHigh          = fread(fid, 1,'float32');
    signal.description(g).ADCLow            = fread(fid, 1,'float32');
    signal.description(g).ADCHigh           = fread(fid, 1,'float32');
    signal.description(g).IndexSignalList   = fread(fid, 1,'int16');
    signal.description(g).CacheOffset       = fread(fid, 1,'int16');
    signal.description(g).Reserved2         = fread(fid,60,'uchar');
    
    % conversion of char values (to right format)
    signal.description(g).SignalName = char(signal.description(g).SignalName(2:signal.description(g).SignalName(1)+1));
    signal.description(g).UnitName   = char(signal.description(g).UnitName(2:signal.description(g).UnitName(1)+1));
  end
  
  % return the header details
  out = signal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read data blocks
NB = signal.header.NumberSampleBlocks;
SD = signal.header.SizeSignalDataBlock;
NS = signal.header.NumberOfSignals;
NS_32bit = NS/2;

if needdat
  
  for g=begblock:endblock
    
    % jump to right position in file
    if signal.header.VersionNumber == 203
      pos = 217 + NS*136 + (g-1) *(86+SD);
    else
      pos = 218 + NS*136 + (g-1) *(86+SD);
    end
    fseek(fid,pos,'bof');
    
    signal.block(g).PI = fread(fid,1,'int32'); %period index
    fread(fid,4,'uchar'); %reserved for extension of previous field to 8 bytes
    signal.block(g).BT = fread(fid,14/2,'int16'); %dostime
    fread(fid,64,'uchar'); %reserved
    data = single(fread(fid,SD/4,'float32'));
    
    % Convert data to 32bit values.
    % In case also 16bit values have to be measured, these values are typecasted below:
    % data = fread(fid,SD/2,'int16'); %read data
    % data = int16(data);
    % data = typecast(data,'int32');
    % signal.block(g).DATA = data;
    signal.data{g} = reshape(data,NS_32bit,SD/(NS_32bit*4));
  end % for
  
  % Converting data to a usable format
  dat = double(cell2mat(signal.data(begblock:endblock)));
  
  for g = 1:NS_32bit  % represent data in [uV]
    dat(g,:) = (dat(g,:) - signal.description(g*2).ADCLow)./(signal.description(g*2).ADCHigh - signal.description(g*2).ADCLow)  .* (signal.description(g*2).UnitHigh - signal.description(g*2).UnitLow) + signal.description(g*2).UnitLow ; %correction for uV
  end
  
  % return the data
  out = dat;
end

fclose(fid);
