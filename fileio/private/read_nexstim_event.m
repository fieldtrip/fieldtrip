function [event] = read_nexstim_event(filename)

% Use as
%   [event] = read_nexstim_event(filename)

% Written by Vladimir Litvak based on the function nxeGetTriggers
% provided by Nexstim
%
% Copyright (C) 2007, Vladimir Litvak
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

