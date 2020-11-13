function [event] = read_nexstim_event(filename)

% Use as
%   [event] = read_nexstim_event(filename)

% Written by Vladimir Litvak based on the function nxeGetTriggers
% provided by Nexstim
%
% Copyright (C) 2007, Vladimir Litvak
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% trigLine - either 1(GATE), 2(TRIG1) or 3(TRIG2)
% trigEdge - either 'rising' or 'falling'

fid=fopen_or_error(filename,'r','l');

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

