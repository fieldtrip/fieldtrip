function [avr] = read_besa_avr(filename)

% READ_BESA_AVR reads average EEG data in BESA format
%
% Use as
%   [avr] = read_besa_avr(filename)
%
% This will return a structure with the header information in
%   avr.npnt
%   avr.tsb
%   avr.di
%   avr.sb
%   avr.sc
%   avr.Nchan   (optional)
%   avr.label   (optional)
% and the ERP data is contained in the Nchan X Nsamples matrix
%   avr.data

% Copyright (C) 2003-2006, Robert Oostenveld
%
% $Log: read_besa_avr.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.4  2008/03/25 10:57:34  roboos
% get channel names from ela file if present
%
% Revision 1.3  2006/10/05 08:53:22  roboos
% added support for another header extension (for Vladimir)
%
% Revision 1.2  2005/04/25 11:10:27  roboos
% changed output labels into cell-array
% added support for reading channel labels from accompanying elp file
%
% Revision 1.1  2005/03/31 07:09:37  roboos
% old implementation, but new implementation of support for avr files that contain channel names
%

fid = fopen(filename, 'rt');

% the first line contains header information
headstr = fgetl(fid);
ok = 0;

if ~ok
  try
    buf = sscanf(headstr, 'Npts= %d TSB= %f DI= %f  SB= %f SC= %f Nchan= %f SegmentName= %s\n');
    avr.npnt  = buf(1);
    avr.tsb   = buf(2);
    avr.di    = buf(3);
    avr.sb    = buf(4);
    avr.sc    = buf(5);
    avr.Nchan = buf(6);
    avr.SegmentName = buf(7);
    ok = 1;
  catch
    ok = 0;
  end
end

if ~ok
  try
    buf = fscanf(headstr, 'Npts= %d TSB= %f DI= %f  SB= %f SC= %f Nchan= %f\n');
    avr.npnt  = buf(1);
    avr.tsb   = buf(2);
    avr.di    = buf(3);
    avr.sb    = buf(4);
    avr.sc    = buf(5);
    avr.Nchan = buf(6);
    ok = 1;
  catch
    ok = 0;
  end
end

if ~ok
  try
    buf = sscanf(headstr, 'Npts= %d TSB= %f DI= %f  SB= %f SC= %f\n');
    avr.npnt  = buf(1);
    avr.tsb   = buf(2);
    avr.di    = buf(3);
    avr.sb    = buf(4);
    avr.sc    = buf(5);
    ok = 1;
  catch
    ok = 0;
  end
end

if ~ok
  error('Could not interpret the header information.');
end

% rewind to the beginning of the file, skip the header line
fseek(fid, 0, 'bof');
fgetl(fid);

% the second line may contain channel names
chanstr = fgetl(fid);
chanstr = deblank(fliplr(deblank(fliplr(chanstr))));
if (chanstr(1)>='A' && chanstr(1)<='Z') || (chanstr(1)>='a' && chanstr(1)<='z')
  haschan   = 1;
  avr.label = str2cell(strrep(deblank(chanstr), '''', ''))';
else
  [root, name] = fileparts(filename);
  haschan = 0;
  elpfile = fullfile(root, [name '.elp']);
  elafile = fullfile(root, [name '.ela']);
  if exist(elpfile, 'file')
    % read the channel names from the accompanying ELP file
    lbl = importdata(elpfile);
    avr.label = strrep(lbl.textdata(:,2) ,'''', '');
  elseif exist(elafile, 'file')
    % read the channel names from the accompanying ELA file
    lbl = importdata(elafile);
    lbl = strrep(lbl ,'MEG ', ''); % remove the channel type
    lbl = strrep(lbl ,'EEG ', ''); % remove the channel type
    avr.label = lbl;
  else
    warning('Could not create channels labels.');
  end
end

% seek to the beginning of the data
fseek(fid, 0, 'bof');
fgetl(fid);   % skip the header line
if haschan
  fgetl(fid); % skip the channel name line
end

buf = fscanf(fid, '%f');
nchan = length(buf)/avr.npnt;
avr.data = reshape(buf, avr.npnt, nchan)';
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to cut a string into pieces at the spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = str2cell(s)
c = {};
[t, r] = strtok(s, ' ');
while ~isempty(t)
  c{end+1} = t;
  [t, r] = strtok(r, ' ');
end
