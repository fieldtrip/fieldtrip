function avr = readBESAavr(filename)
% READBESAAVR reads sensor level data from an AVR-file. 
%
% Parameters:
%     [filename]
%         In the case that the current folder is not the folder containing 
%         the file it should be the full path including the name of the 
%         elp file else only the name of the file should be specified. 
% 
% Return:
%     [avr] 
%         A Matlab structure containing the data and the corresponding
%         parameters stored in the AVR-file.
% 
% Copyright (C) 2013, BESA GmbH
%
% File name: readBESAavr.m
%
% Author: Todor Jordanov
% Created: 2013-07-30

fp = fopen(filename, 'r');

% MATLAB reserves file identifiers 0, 1, and 2 for standard input,  
% standard output (the screen), and standard error, respectively. When 
% fopen successfully opens a file, it returns a file identifier greater 
% than or equal to 3.
if(fp >= 3)
    
    % Get the first line of the file. It looks something like that:
    % Npts= 512   TSB= -400.000 DI= 3.125000 SB= 1.000 SC= 200.0 ...
    FirstLine = fgetl(fp);
    % First of all check if the parameter SegmentName exists.
    if(~isempty(strfind(FirstLine, 'SegmentName')))
        
        tmp = regexp(FirstLine, 'SegmentName= ', 'split');
        avr.SegmentName = tmp{2};
        FirstLine = tmp{1};
        
    end
    
    % Parse header information of .avr file.
    [DataParams] = sscanf(FirstLine, ...
        'Npts= %i TSB= %f DI= %f SB= %f SC= %f Nchan= %i');
    avr.Npts = DataParams(1);
    avr.TSB = DataParams(2);    
    avr.DI = DataParams(3);
    avr.Time = avr.TSB:avr.DI:avr.TSB+avr.DI*(avr.Npts-1);
    
    if(size(DataParams, 2) > 5)
        
        avr.Nchan = DataParams(6);
        
    end
    
    avr.ChannelLabels = [];
    
    % The second line could contain the channel labels but it also could be
    % that it contains data.
    SecondLine = fgetl(fp);
    
    % Check if the second line contains labels or data values
    Characters = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' ...
        'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};
    NumChars = size(Characters, 2);
    iCounter = 1;
    bContainsCharacter = 0;
    for ch = 1:NumChars
        
        % If the second line contains characters then these are the channel
        % labels.
        if(~isempty(regexpi(SecondLine, Characters{ch}, 'match')))
            
            avr.ChannelLabels = SecondLine;
            bContainsCharacter = 1;
            break;
        
        end
        
    end
    
    % If the second line does not contain characters then these are no
    % labels but data instead.
    if(bContainsCharacter == 0)
        
        avr.Data(iCounter, :) = sscanf(SecondLine, '%f', [avr.Npts,1]);
        iCounter = 2;
        
    end

    while(true)
        
        CurrentLine = fgetl(fp);
        % Check if end of file.
        if(~ischar(CurrentLine))
            
            break;
            
        end
        
        avr.Data(iCounter, :) = sscanf(CurrentLine, '%f', [avr.Npts,1]);
        iCounter = iCounter + 1;
        
    end
    
    avr.Nchan=size(avr.Data, 1);
    
    fclose(fp);
    
else
    
    avr = [];
    disp('Error! Invalid file identifier.')
    
end