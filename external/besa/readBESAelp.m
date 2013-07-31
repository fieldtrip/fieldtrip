function Coords = readBESAelp(filename)
% READBESAELP reads the channel position from an ELP-file. 
%
% Parameters:
%     [filename]
%         In the case that the current folder is not the folder containing 
%         the file it should be the full path including the name of the 
%         elp file else only the name of the file should be specified. 
% 
% Return:
%     [Coords] 
%         A matrix containing the polar coordinates of the electrodes. The
%         size of the matrix is [NumChannels x 2].
% 
% Copyright (C) 2013, BESA GmbH
%
% File name: readBESAelp.m
%
% Author: Todor Jordanov
% Created: 2013-07-30

fp = fopen(filename, 'r');

if (fp)

    % Read the first line of the file. It should look like this:
    % [EEG] [Fp1] -92 -72 [1]
    % the parameters in [] are not always available.
    FirstLine = fgetl(fp);
    
    ChannelCounter = 1;
    
    tmp = regexp(FirstLine, ' ', 'split');
    
    % Get the number of the columns.
    NumParams = size(tmp, 2);
    
    IndexPhi = 0;
    IndexTheta = 0;
    
    % Check how many parameters are available and set the indices of the
    % coordinate values correspondingly.
    switch NumParams
        case 2
            IndexPhi = 1;
            IndexTheta = 2;
        case 3
        	IndexPhi = 2;
            IndexTheta = 3;
        case 4
        	IndexPhi = 3;
            IndexTheta = 4;
        case 5
        	IndexPhi = 3;
            IndexTheta = 4;
        otherwise
        	disp('Wrong number of parameters!!!')
    end

    if(IndexPhi ~= 0 && IndexTheta ~=0)
        
        Coords(1, 1) = str2double(tmp(IndexPhi));
        Coords(1, 2) = str2double(tmp(IndexTheta));
        
        while(true)
            
            CurrentLine = fgetl(fp);
            % Check if end of file.
            if(~ischar(CurrentLine))

                break;

            end
            
            ChannelCounter = ChannelCounter + 1;
            tmp = regexp(CurrentLine, ' ', 'split');
            Coords(ChannelCounter, 1) = str2double(tmp(IndexPhi));
            Coords(ChannelCounter, 2) = str2double(tmp(IndexTheta));
            
        end
        
    end
    
end