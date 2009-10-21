function [rej] = read_brainvision_marker(fn);

% READ_BRAINVISION_MARKER reads rejection marks from a BrainVision file
%
% Use as
%   [rej] = read_brainvision_marker(filename)
%
% This function returns a Nx2 matrix with the begin and end latency
% of N rejection marks. The latency is in miliseconds.

% Copyright (C) 2004, Robert Oostenveld & Doug Davidson
%
% $Log: read_brainvision_marker.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.2  2004/09/24 15:57:20  roboos
% implemented suggested change by Doug Davidson: get sampling rate from file instead of having it hard coded
%
% Revision 1.1  2004/08/26 15:53:18  roboos
% new implementation by Doug Davidson, based upon read_eep_rej
%

rej = [];

[x y] = textread(fn, '%s%s', 1, 'delimiter', ',:' );

p   = mat2str(cell2mat(y));
pos = findstr('Hz', p);

samplingrate = str2num(p(1:pos-1));

s = 1./samplingrate;

[col1 col2] = textread(fn, '%*s%*s%n%n%*[^\n]',...
    'delimiter', ',',...
    'headerlines', 2 );

% Start and end points in msec
startp = (col1*s)*1000;
endp = ((col2*s)+(col1*s))*1000;

rej = [startp endp];
clear x p pos samplingrate;
