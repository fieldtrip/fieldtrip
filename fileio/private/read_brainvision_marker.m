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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
