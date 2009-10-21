function [swf] = read_besa_swf(filename);

% READ_BESA_SWF
%
% Use as
%   [swf] = read_besa_swf(filename)
%
% This will return a structure with the header information in
%   swf.label     cell array with labels
%   swf.data      data matrix, Nchan X Npnts
%   swf.npnt
%   swf.tsb
%   swf.di
%   swf.sb

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: read_besa_swf.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.3  2006/06/22 15:03:07  roboos
% fixed bug, data was transposed in case of row-formatted file
%
% Revision 1.2  2006/03/20 08:39:13  roboos
% 2 small bug fixes, thanks to Vladimir
%
% Revision 1.1  2006/03/16 17:31:42  roboos
% new implementation, supports both rows and columns
%

fid = fopen(filename);
line = fgetl(fid);
head = sscanf(line,'Npts= %f TSB= %f DI= %f SB= %f'); % new BESA versions include more information in the .swf file header; skip that
Npts = head(1);
TSB  = head(2);
DI   = head(3);
SB   = head(4);

offset = ftell(fid);

% try reading a label and a number
dum = textscan(fid, '%s %f', 1);
fseek(fid, offset, 'bof');

if  isempty(dum{2})
  % first line contains all labels, each subsequent line is one timepoint for all channels
  line = fgetl(fid);
  line(line==' ') = [];
  label = tokenize(line, ':')';
  if isempty(label{end})
    label = label(1:(end-1));  % the last one is empty
  end
  Nsrc = length(label);
  dat = fscanf(fid, '%f', inf);
  dat = reshape(dat, Nsrc, Npts);
else
  % each line starts with a single label, followed by all timepoints for that channel
  dat = [];
  label = {};
  count = 0;
  while ~feof(fid)
    count = count+1;
    % read one line at a time
    line = fgetl(fid);
    sel = find(line==':', 1);
    label{count} = deblank(line(1:(sel-1)));
    buf = sscanf(line((sel+1):end), '%f');
    dat(count,:) = buf(:);
  end
end

fclose(fid);

% assign the output, should be similar as READ_BESA_AVR
swf.npnt  = Npts;
swf.tsb   = TSB;
swf.di    = DI;
swf.sb    = SB;
swf.label = label;
swf.data  = dat;

% FIXME, here combine the channels of the regional sources


