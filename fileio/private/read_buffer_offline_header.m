function [hdr, nameFlag] = read_buffer_v1_header(headerfile)
% function [hdr, nameFlag] = read_buffer_v1_header(headerfile)
%
% This function reads a FCDC buffer header from a binary file.
%
% On return, nameFlag has one of the following values:
%   0 = No labels were generated (fMRI etc.)
%   1 = Fake labels were generated
%   2 = Got channel labels from chunk information

% (C) 2010 S. Klanke

type = {
  'char'
  'uint8'
  'uint16'
  'uint32'
  'uint64'
  'int8'
  'int16'
  'int32'
  'int64'
  'single'
  'double'
};

wordsize = {
  1 % 'char'
  1 % 'uint8'
  2 % 'uint16'
  4 % 'uint32'
  8 % 'uint64'
  1 % 'int8'
  2 % 'int16'
  4 % 'int32'
  8 % 'int64'
  4 % 'single'
  8 % 'double'
};

hdr = [];
hdr.nSamples    = 0;
%hdr.nEvents     = 0;
hdr.nSamplesPre = 0;
hdr.nTrials     = 1;
hdr.orig =[];

nameFlag = 0;

FA = fopen([headerfile '.txt'], 'r');
while ~feof(FA);
  s = fgetl(FA);
  indEq = find(s=='=');
  if numel(indEq)==1
    key = s(1:(indEq-1));
    val = s((indEq+1):end);
    switch key
      case 'fSample'
        hdr.Fs = str2num(val);
      case 'nSamples'
        hdr.nSamples = str2num(val);
      case 'nEvents'
%        hdr.nEvents = str2num(val);
      case 'nChans'
        hdr.nChans = str2num(val);
        hdr.label = cell(hdr.nChans,1);
      case 'dataType'
        if strcmp(val,'float32')
           hdr.orig.data_type = 'single';
        elseif strcmp(val,'float64')
           hdr.orig.data_type = 'double';
        else
          hdr.orig.data_type = val;
        end
      case 'version'
        hdr.orig.version = str2num(val);
      case 'endian'
        if strcmp(val,'big')
          hdr.orig.endianness='ieee-be';
        elseif strcmp(val,'little')
          hdr.orig.endianness='ieee-le';
        else
          error('Invalid value for key ''endian''.');
        end
    end
    continue;
  end
  indCol = find(s==':');
  if numel(indCol)>=1
    ind = str2num(s(1:(indEq-1)));
    name = s((indEq+1):end);
    if ~isempty(ind) && ind>=1 && ind<=hdr.nChans
      hdr.label{ind} = name;
      nameFlag = 2;
    else
      error('Invalid channel name definition - broken header file?');
    end
  end
end

fclose(FA);

F = fopen(headerfile,'rb',hdr.orig.endianness);
bin=[];
bin.nChans    = double(fread(F,1,'uint32'));
bin.nSamples  = double(fread(F,1,'uint32'));
bin.nEvents   = double(fread(F,1,'uint32'));
bin.Fs        = double(fread(F,1,'single'));
bin.data_type = fread(F,1,'uint32');
bin.bufsize   = fread(F,1,'uint32');

if hdr.nChans ~= bin.nChans
   error('Number of channels in binary header does not match ASCII definition');
end
if hdr.nSamples ~= bin.nSamples
   warning('Number of samples in binary header does not match ASCII definition');
end
%if hdr.nEvents ~= bin.nEvents
%   error('Number of samples in binary header does not match ASCII definition');
%end
if strcmp(hdr.orig.data_type, type{bin.data_type+1})
   hdr.orig.wordsize = wordsize{bin.data_type+1};
else
   error('Data type in binary header does not match ASCII definition');
end

if hdr.nSamples == 0
  datafile = [headerfile(1:end-6) 'samples'];
  FD = fopen(datafile, 'rb');
  fseek(FD,0,'eof');
  fileSize = ftell(FD);
  fclose(FD);
  
  n = fileSize / hdr.orig.wordsize / hdr.nChans;
  hdr.nSamples = floor(n);
  if hdr.nSamples ~= n
    warning('Size of ''samples'' is not a multiple of the size of one sample');
  end
end

% TODO: add chunk handling
while ~feof(F)
	typeAndSize = fread(F, 2, 'uint32');
	if numel(typeAndSize) < 2
		break
	end
	switch typeAndSize(1)
		case 1 % channel names
			% already dealt with, TODO: maybe check consistency with ASCII stuff
			dummy = fread(F, typeAndSize(2), 'uint8=>uint8');
		case 2 % FT_CHUNK_CHANNEL_FLAGS 
			% FIXME: ignored for now
			dummy = fread(F, typeAndSize(2), 'uint8=>uint8');
		case 3 % FT_CHUNK_RESOLUTIONS 
			if typeAndSize(2) == 8*hdr.nChans
				hdr.resolutions = fread(F, [hdr.nChans, 1], 'double');
			else
				warning('Invalid size of RESOLUTIONS chunk - skipping');
				dummy = fread(F, typeAndSize(2), 'uint8=>uint8');
			end
		case 4 % FT_CHUNK_ASCII_KEYVAL
			dummy = fread(F, typeAndSize(2), 'uint8');
			% FIXME: ignore for now
		case 5 % FT_CHUNK_NIFTI1
			hdr.orig.nifti_1 = fread(F, [1, typeAndSize(2)], 'uint8=>uint8');
			if typeAndSize(2) == 348
				hdr.nifti_1 = decode_nifti1(hdr.orig.nifti_1);
			else
				warning('Invalid size of NIFTI_1 chunk - skipping');
			end
		case 6 % FT_CHUNK_SIEMENS_AP = 6
			hdr.orig.siemensap = fread(F, typeAndSize(2), 'uint8=>uint8');
			if exist('sap2matlab')==3
				hdr.siemensap = sap2matlab(hdr.orig.siemensap);
			end
		case 7 % FT_CHUNK_CTF_RES4 = 7
			hdr.orig.ctf_res4 = fread(F, typeAndSize(2), 'uint8=>uint8');
			tmp_name = tempname;
			FT = fopen(tmp_name, 'wb');
			fwrite(FT, orig.ctf_res4, 'uint8');
			fclose(FT);
			R4F = read_ctf_res4(tmp_name);
			delete(tmp_name);
			% copy over the labels
			hdr.label = R4F.label;
			% copy over the 'original' header
			hdr.orig = R4F;
	end
end

fclose(F);

if nameFlag < 2 && hdr.nChans < 2000
  nameFlag = 1; % fake labels generated - these are unique
  for i=1:hdr.nChans
    hdr.label{i} = sprintf('%d', i);
  end
end

