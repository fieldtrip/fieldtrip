function dat = read_buffer_offline_data(datafile, hdr, range)
% function dat = read_buffer_offline_data(datafile, header, range)
%
% This function reads FCDC buffer-type data from a binary file.

% (C) 2010 S. Klanke

type2type = [hdr.orig.data_type '=>' hdr.orig.data_type];
nStart = range(1)-1; % 0-based offset in samples
bStart = nStart * hdr.orig.wordsize * hdr.nChans;
nRead = range(2)+1-range(1);
bRead = nRead * hdr.orig.wordsize * hdr.nChans;

k = 0;
sizeSoFar = 0;
% start with name_k = datafile (= '.../samples' )
name_k = datafile;

while true
  F = fopen_or_error(datafile,'rb',hdr.orig.endianness);
  fseek(F, 0, 'eof');
  size_k = ftell(F);
   
  if nStart >= sizeSoFar && nStart < sizeSoFar + size_k
    fseek(F, bStart - sizeSoFar, 'bof');
    if bStart + bRead <= sizeSoFar + size_k
      % desired region is completely contained in this file
	  dat = fread(F,[hdr.nChans nRead], type2type);
	  bReadRemain = 0;
	else
	  % desired region starts in this file, but then goes on
	  bRead_k = sizeSoFar + size_k - bStart;
	  nRead_k = bRead_k / (hdr.orig.wordsize * hdr.nChans);
	  
	  dat = fread(F,[hdr.nChans nRead_k], type2type);
	  bReadRemain = nRead - nRead_k;
	end
  else
    fseek(F, 0, 'bof');
	if bReadRemain <= size_k
	  % this is the last file we need
	  bRead_k = bReadRemain;
	  bReadRemain = 0;
	else 
	  bRead_k = size_k;
	  bReadRemain = bReadRemain - size_k;
	end
	nRead_k = bRead_k / (hdr.orig.wordsize * hdr.nChans);
	dat = [dat fread(F,[hdr.nChans nRead_k], type2type)];
  end	  
  fclose(F);		 
  if bReadRemain == 0
     break;
  end
  % update name_k for next iteration
  k=k+1;
  name_k = sprintf('%s%i', datafile, k);
end
