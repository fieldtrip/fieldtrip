function [event] = read_nexstim_event(filename)

% Use as
%   [event] = read_nexstim_event(filename)

% Written by Vladimir Litvak based on the function nxeGetTriggers
% provided by Nexstim
%
% Copyright (C) 2007, Vladimir Litvak
%
% $Log: read_nexstim_event.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.2  2007/12/17 13:03:52  roboos
% Vladimir found and fixed some bugs pertaining to the nexstim_nxe format
%
% Revision 1.1  2007/12/17 08:24:28  roboos
% added support for nexstim_nxe, thanks to Vladimir
% the low-level code has not been tested by myself
%

% trigLine - either 1(GATE), 2(TRIG1) or 3(TRIG2)
% trigEdge - either 'rising' or 'falling'

fid=fopen(filename,'r','l');

numChannels = 64;
blockSamples = 14500;
trigThreshold = 2000;
trigChannels = [1 2 3]; % Is fixed for the present Nexstim format.

fseek(fid,0,'eof');
numBytes = ftell(fid);
numSamples = (numBytes/2)/numChannels;
numBlocks = ceil(numSamples/blockSamples);

fseek(fid,0,'bof');
trigPos=[];

event=[];
for i = 1:numBlocks
    blockPos=(i-1)*blockSamples;
    data = fread(fid,[numChannels blockSamples+1],'int16');

    for trigLine=1:length(trigChannels);
        trigPosRising = find(diff(data(trigChannels(trigLine),:))>trigThreshold)+blockPos;
        trigPosFalling = find(diff(data(trigChannels(trigLine),:))<-trigThreshold)+blockPos;

        if ~isempty(trigPosRising)
            for i=1:length(trigPosRising)
                event(end+1).type     = 'rising';
                event(end  ).sample   =  trigPosRising(i);
                event(end  ).value    = trigLine;
                event(end  ).offset   = [];
                event(end  ).duration = [];
            end
        end

        if ~isempty(trigPosFalling)
            for i=1:length(trigPosFalling)
                event(end+1).type     = 'falling';
                event(end  ).sample   =  trigPosFalling(i);
                event(end  ).value    = trigLine;
                event(end  ).offset   = [];
                event(end  ).duration = [];
            end
        end

    end
    fseek(fid,-(2*numChannels),'cof');
end

fclose(fid);

