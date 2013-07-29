function swf = readBESAswf(filename)

% readBESAswf read all information from an *.swf file
%
% Note that the *.swf file must contain waveforms in rows and not in
% columns.
%
% Use as
%   swf = readBESAswf(filename)
%
% The output is a structure containing the following fields:
%   Npts: Number of time points
%   Time: array of sampled time instants
%   wavename: Labels of the different sources
%   data: matrix of data points [Number of sources x Npts]
%
% Modified April 26, 2006 Robert Oostenveld
% Modified June 27, 2006 Karsten Hoechstetter: Function reads now both
%      column and row swf formats

if isempty(findstr(filename,'.'))
    filename = [filename,'.swf'];
end
fp = fopen(filename);

% Read the file
if (fp)
    % Read header of .swf file
    headline = fscanf(fp,'Npts= %f TSB= %f DI= %f SB= %f');
    % Read required information.
    %swf.Npts = fscanf(fp,'Npts= %i');
    %swf.TSB = fscanf(fp,' TSB= %f');    
    %swf.DI = fscanf(fp,'DI= %f');
    
    %try 
    %    swf.SB = fscanf(fp,'SB= %f'); 
    %catch
    %end
    
    % New BESA versions include more information in the .swf file header; skip that
    i=1;
    while i<1000
        a = fscanf(fp,'%c',1);
        if strcmp(a,sprintf('\n'))
        i=1000;
        end
        i=i+1;
    end
    swf.Npts=headline(1);
    TSB=headline(2);
    DI=headline(3);
    swf.Time = [TSB:DI:TSB+(swf.Npts-1)*DI];

    % Read first line after header to decide whether swf's are in rows or
    % columns
    row=1;
    SecondLine = fgets(fp);
    a=strfind(SecondLine,': ');
    if a(end) >= length(SecondLine)-10      %swf's are in columns (because there is a ":" within the last 11 characters of the second line)
        row=0;
    end

    % Read data and labels
    if row == 1             % swf's in rows
        Name = sscanf(SecondLine,'%s',1);
        swf.data(1,:) = sscanf(SecondLine(length(Name)+1:end),'%f',[1 headline(1)]);
        swf.waveName(1) = cellstr(Name(1:length(Name)-1));
        for i=2:1000
            try             % check if there is another waveform
                Name = fscanf(fp,'%s',1);            
                swf.data(i,:) = fscanf(fp,'%f',[1,headline(1)]); 
                swf.waveName(i) = cellstr(Name(1:length(Name)-1));
            catch           % stop if end of file is reached
                break  
            end
        end
    else                    % swf's in columns
        temp1 = SecondLine(1:a(1)-1);              
        temp2 = deblank(temp1(end:-1:1));
        swf.waveName(1) = cellstr(temp2(end:-1:1));
        for i=2:length(a)
            temp1 = SecondLine(a(i-1)+2:a(i)-1);
            temp2 = deblank(temp1(end:-1:1));
            swf.waveName(i) = cellstr(temp2(end:-1:1));
        end
        swf.data = fscanf(fp,'%f',[length(a),headline(1)]);
        i=length(a)+1;
    end
end

fclose(fp);