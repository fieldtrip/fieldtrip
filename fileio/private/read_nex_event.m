function [event] = read_nex_event(filename)

% READ_NEX_EVENT for Plexon *.nex file
%
% Use as
%   [event] = read_nex_event(filename)
%
% The sample numbers returned in event.sample correspond with the
% timestamps, correcting for the difference in sampling frequency in the
% continuous LFP channels and the system sampling frequency. Assuming 40kHz
% sampling frequency for the system and 1kHz for the LFP channels, it is
%   event.sample = timestamp / (40000/1000);
%
% See also READ_NEX_HEADER, READ_NEX_DATA

% Copyright (C) 2005-2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

hdr = read_nex_header(filename);
adindx = find(cell2mat({hdr.varheader.typ})==5);
smpfrq = hdr.varheader(adindx(1)).wfrequency;

% find the channel with the strobed trigger
mrkvarnum = find([hdr.varheader.typ] == 6);

fid=fopen(filename,'r','ieee-le');
status = fseek(fid,hdr.varheader(mrkvarnum).offset,'bof');

% read the time of the triggers
dum = fread(fid,hdr.varheader(mrkvarnum).cnt,'int32');
dum = dum ./(hdr.filheader.frequency./smpfrq);
mrk.tim = round(dum);

% read the value of the triggers
status = fseek(fid,64,'cof');
dum = fread(fid,[hdr.varheader(mrkvarnum).mrklen,hdr.varheader(mrkvarnum).cnt],'uchar');
mrk.val = str2num(char(dum(1:5,:)'));

status = fclose(fid);

% translate into an FCDC event structure
Nevent = length(mrk.tim);
event = struct('sample', num2cell(mrk.tim), 'value', num2cell(mrk.val));
for i=1:Nevent
  % the code above with the struct and num2cell is much faster
  % event(i).sample   = mrk.tim(i);
  % event(i).value    = mrk.val(i);
  event(i).type     = hdr.varheader(mrkvarnum).nam;
  event(i).duration = 1;
  event(i).offset   = 0;
end

