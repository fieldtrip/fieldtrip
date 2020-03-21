function mat2ades(dat, fileName, Fs, chanlabel, chanunit, dattype, datunit)

% write in the current folder ADES and DAT files from matrix in MATLAB workspace
% data = matrix of data (nbchannel * time points) - the data have to be in microVolt
% fileName = string of the output files without extension ; the ades and dat files will have the same name
% FS = sampling rate
% labels = cell-array with channel labels
% labelType : 'EEG' or 'MEG'
%
% Data are stored in a binary file which name is exactly the same than the header file except the extension: .dat
% The samples are stored as float, 4 bytes per sample, little endian. The channels are multiplexed.
%
% Sophie Chen - January 2014
% Modified by Robert Oostenveld - February 2019

%% generate the ADES file
adesFile = [fileName '.ades'];

fid = fopen_or_error(adesFile, 'wt');

fprintf(fid, '#ADES header file\r\n');
fprintf(fid, 'samplingRate = %d\r\n', Fs);
fprintf(fid, 'numberOfSamples = %d\r\n', size(dat,2));

for lab = 1:length(dattype)
  fprintf(fid, 'Unit = %s, %s\r\n', dattype{lab}, datunit{lab});
end

for lab = 1:length(chanlabel)
  fprintf(fid, '%s = %s\r\n', chanlabel{lab}, chanunit{lab});
end

fclose(fid);

%% generate the DAT file

datFile = [fileName '.dat'];

fad = fopen_or_error(datFile, 'wb');
fwrite(fad, dat, 'float32', 'l');
fclose(fad);
end