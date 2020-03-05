function write_gdf(filename, hdr, data)

% WRITE_GDF(filename, header, data)
%
% Writes a GDF file from the given header (only label, Fs, nChans are of interest)
% and the data (unmodified). Digital and physical limits are derived from the data
% via min and max operators. The GDF file will contain N records of 1 sample each,
% where N is the number of columns in 'data'.
  
% Copyright (C) 2010, Stefan Klanke
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

[nChans,N] = size(data);
if hdr.nChans ~= nChans
  error 'Data dimension does not match header information';
end
fSample = real(hdr.Fs(1)); % make sure this is a scalar + real

labels = char(zeros(nChans, 16));
% labels
for n=1:nChans
  ln = length(hdr.label{n});
  if ln > 16
    fprintf(1, 'Warning: truncating label %s to %s\n', hdr.label{n}, hdr.label{n}(1:16));
    ln = 16;
  end
  labels(n,1:ln) = hdr.label{n}(1:ln);
end

if ~isreal(data)
  error 'Cannot write complex-valued data.';
end

switch class(data)
  case 'int8'
    gdfType = 1;
  case 'uint8'
    gdfType = 2;
  case 'int16'
    gdfType = 3;
  case 'uint16'
    gdfType = 4;
  case 'int32'
    gdfType = 5;
  case 'uint32'
    gdfType = 6;
  case 'int64'
    gdfType = 7;
  case 'uint64'
    gdfType = 8;
  case 'single'
    gdfType = 16;
  case 'double'
    gdfType = 17;
  otherwise
    error 'Unrecognized data type';
end

% TODO: check if this makes sense
% will be terrible for appending data...
digMin = double(min(data,[],2));
digMax = double(max(data,[],2));

% adjust the the digital min/max a bit, otherwise the biosig reading code
% will return NaNs for the most extreme values
digMin = digMin - 1e5.*eps(digMin);
digMax = digMax + 1e5.*eps(digMax);

physMin = digMin;
physMax = digMax;

%struct GDF_Header {
%	char version[8];
%	char patient[66];
% uint8_t reserved[10];
%	uint8_t smokingAlcoholDrugMedication;
%	uint8_t weight;
%	uint8_t height;
%	uint8_t genderHandedVisualHeart;
%	char recordingId[64];
%	uint32_t recordingLocation[4];
%	uint32_t recordingTime[2];
%	uint32_t birthday[2];
%	uint16_t headerLengthInBlocks;
%	uint8_t patientClass[6];
%	uint64_t equipmentId;
%	uint8_t reserved2[6];
%	uint16_t headsize[3];
%	float posRefElec[3];
%	float posGroundElec[3];
%	int64_t numDataRecords;
%	uint32_t durDataRecord[2];
%	uint16_t numChannels;
%	uint16_t reserved3;
  
fid = fopen_or_error(filename, 'wb', 'ieee-le');
% first write fixed part
fprintf(fid, 'GDF 2.20'); %version
fwrite(fid, zeros(1,66), 'int8'); % patient
fwrite(fid, zeros(1,10+4), 'uint8'); % reserved + smoking/weight/height/gender
fwrite(fid, zeros(1,64), 'int8');    % recording ID
fwrite(fid, zeros(1,4), 'uint32');   % recording location
secSince2000 = etime(clock, [2000 1 1 0 0 0]);
fwrite(fid, secSince2000, 'uint32'); % recording time in seconds since 2000/1/1
fwrite(fid, 0, 'int32'); % fractional part ignored

fwrite(fid, zeros(1,2), 'uint32'); % birthday
hdrLengthInBlocks = nChans + 1; % 256 bytes = 1 block per channel + 256 bytes for fixed header
fwrite(fid, hdrLengthInBlocks, 'uint16');
fwrite(fid, zeros(1,6), 'uint8'); % patient class
fwrite(fid, 0, 'uint64'); % equipment ID
fwrite(fid, zeros(1,6), 'uint8'); % reserved2
fwrite(fid, zeros(1,3), 'uint16'); % headsize
fwrite(fid, zeros(1,6), 'single'); % posRefElec + posGroundElec

fwrite(fid, N, 'int64'); % numDataRecords = #samples
[num,den] = rat(1.0/fSample); % sample rate as rational number
fwrite(fid, num, 'uint32');   % nominator
fwrite(fid, den, 'uint32');   % denominator
fwrite(fid, nChans, 'uint16'); % number of channel
fwrite(fid, zeros(1,1), 'uint16'); % reserved3

% now variable part of header, 256 bytes per channel
fwrite(fid, labels', 'char*1');
fwrite(fid, zeros(80,nChans), 'uint8'); % Types
fwrite(fid, zeros(6,nChans), 'uint8'); % PhysDim
fwrite(fid, zeros(1,nChans), 'uint16'); % PhysDimCode

fwrite(fid, physMin, 'double'); % physical minimum
fwrite(fid, physMax, 'double'); % physical maximum
fwrite(fid, digMin, 'double');  % digital minimum
fwrite(fid, digMax, 'double');  % digital maximum
fwrite(fid, zeros(68,nChans), 'uint8'); % pre-filter information
fwrite(fid, zeros(nChans, 3), 'single'); % lowpass, highpass, notch
fwrite(fid, ones(1, nChans), 'int32'); % samples per record = 1
fwrite(fid, gdfType*ones(1, nChans), 'uint32'); % gdf data type
fwrite(fid, zeros(nChans, 3), 'single'); % sensor position
fwrite(fid, zeros(nChans, 20), 'uint8'); % sensor description

fwrite(fid, data, class(data));
fclose(fid);
