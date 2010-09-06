function dat = read_buffer_offline_data(datafile, hdr, range)
% function dat = read_buffer_offline_data(datafile, header, range)
%
% This function reads FCDC buffer-type data from a binary file.

% (C) 2010 S. Klanke

F = fopen(datafile,'rb',hdr.orig.endianness);

nStart = range(1)-1; % 0-based offset in samples
nStart = nStart * hdr.orig.wordsize * hdr.nChans;
nRead = range(2)+1-range(1);

type2type = [hdr.orig.data_type '=>' hdr.orig.data_type];
fseek(F, nStart, 'bof');
dat = fread(F,[hdr.nChans nRead], type2type);
fclose(F);