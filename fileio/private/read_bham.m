function [dat, lab] = read_bham(filename)

% READ_BHAM reads the EEG data files as recorded by Praamstra in Birmingham
% the datafiles are in a particular ascii format
%
% [dat, lab] = read_bham(filename)

% Copyright (C) 2000, Robert Oostenveld
% 
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fid = fopen(filename, 'rt');

lablen = 6;
line   = fgetl(fid);
numelc = 0;
while ~isempty(line)
  numelc = numelc + 1;
  [t, r] = strtok(line);
  lab(numelc,:) = [blanks(lablen-length(t)), t];
  line = r;
end

buf = fscanf(fid, '%f');
dat = zeros(numelc, length(buf)/numelc);
dat(:) = buf;

