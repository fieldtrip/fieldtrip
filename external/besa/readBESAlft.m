function [dim lf] = readBESAlft(filename, chanid)

% readBESAlft reads the leadfield from the BESA leadfield file format.
% The leadfield matrix and a vector with the leadfield dimensions are
% returned.
%
% Parameters:
%     [filename]
%         In the case that the current folder is not the folder containing 
%         the file it should be the full path including the name of the 
%         lft file else only the name of the file should be specified. 
% 
%     [chanid]
%         Additional parameter specifying particular channel for which
%         leadfield information should be returned.
%
% Return:
%     [dim] 
%         Leadfield dimensions {<number sensors>, <number source space
%         nodes>, <number source directions>}. Please note, in case that
%         the additional parameter is specified, the number of sensors will
%         be set to 1.
% 
%     [lf] 
%         Leadfield matrix with <number sensors> rows and <number source
%         space nodes> x <number source directions> columns. Each row
%         contains first the potentials for sources at all nodes in x-dir,
%         then for sources at all nodes in y-dir, and finally for all
%         nodes in z-dir.
%
% Copyright (C) 2013, BESA GmbH
%
% File name: readBESAlft.m
%
% Author: Benjamin Lanfer/Robert Spangler
% Created: 2013-11-27


% Check additional parameter
if exist('chanid', 'var')
    % check if chanid is scalar
    if isscalar(chanid)
        % get leadfield info from specific channels
    else
        % error; only scalar values allowed for chanid
        error('Parameter chanid is not a scalar value!');
    end;
else
    % get leadfield info all channels
    chanid = [];
end;

% Try to open file.
FileID = fopen(filename, 'rb');
if(FileID < 0)
	fprintf('Error opening file.\n');
	return
end

% Read version number (int32)
[VersionNumber NrReadElements] = fread(FileID, 1, 'int32');
if(NrReadElements ~= 1)
	fprintf('Could not read number of elements.\n');
	return
end

% Check version number
ExpectedVersionNumber = 1;
if(VersionNumber ~= ExpectedVersionNumber)
	fprintf('Wrong version number. Expected: %d, read %d.\n', ExpectedVersionNumber, ...
		VersionNumber);
	return
end

% Read number of sensors (int32)
[NumberSensors NrReadElements] = fread(FileID, 1, 'int32');
if(NrReadElements ~= 1)
	fprintf('Could not read number of sensors.\n');
	return
end

% In case the leadfield for one specific channel has to be returned, check
% if chanid is in correct range.
if ~isempty(chanid)
    if ((chanid > NumberSensors) ||...
        (chanid < 1))
        % error; chanid out of range
        error('Parameter chanid is out of range!');    
    end;
end;

% Read number of source space nodes (int32)
[NumberSourceSpaceNodes NrReadElements] = fread(FileID, 1, 'int32');
if(NrReadElements ~= 1)
	fprintf('Could not read number of source space nodes.\n');
	return
end

% Read number of source directions per node (int32)
[NumberDirections NrReadElements] = fread(FileID, 1, 'int32');
if(NrReadElements ~= 1)
	fprintf('Could not read number of source directions.\n');
	return
end

% Read maximum LF values for each source (source nodes x directions) (float32)
NumberColumns = NumberDirections * NumberSourceSpaceNodes;
[MaxVals NrReadElements] = fread(FileID, NumberColumns, 'float32');
if(NrReadElements ~= NumberColumns)
	fprintf('Could not read maximum leadfield values from file.\n');
	return
end

if ~isempty(chanid)
    % Get leadfield for specific channel only
    lf = [];
    NumberSensors = 1;
    
    % Go to correct position in file where leadfield info for current
    % channel is stored
    NumSkippedBytes = (chanid-1)*NumberSourceSpaceNodes*NumberDirections*2;
    status = fseek(FileID, NumSkippedBytes, 0);
    if status ~= 0
        error('Unable to move to specified position in file.');
    end;

    % Read compactly stored LF values for channel #chanid (int16)
    NumberLFValues = NumberColumns * 1;
    [CompactLFVals NrReadElements] = fread(FileID, NumberLFValues, 'int16');
    if(NrReadElements ~= NumberLFValues)
        fprintf('Could not read leadfield values from file.\n');
        return
    end
    
else
    % Get leadfield for all channels

    % Read compactly stored LF values (int16)
    NumberLFValues = NumberColumns * NumberSensors;
    [CompactLFVals NrReadElements] = fread(FileID, NumberLFValues, 'int16');
    if(NrReadElements ~= NumberLFValues)
        fprintf('Could not read leadfield values from file.\n');
        return
    end
end

% Reshape matrix with compact LF values.
CompactLFVals = reshape(CompactLFVals, NumberColumns, NumberSensors);
CompactLFVals = CompactLFVals';

% Compute factors for converting compactly stored LF values
% to float values.
ConversionFactors = MaxVals / 32767.0;

% Clear superfluous parameters
clear('MaxVals');

% Convert compact LF values to float.
%lf = CompactLFVals * repmat(ConversionFactors, 1, NumberColumns);

% NOTE RS:
% * operator in line %lf = CompactLFVals * repmat(...) is a real
% matrix multiplication, not a scalar multiplication. Is this intended?

% NOTE RS:
% Use alternative conversion of compact LF values to float LF
% values due to the extrem memory consumption of the repmat function:
for iCol=1:NumberColumns
    CompactLFVals(:,iCol) = CompactLFVals(:,iCol)*ConversionFactors(iCol);
end
lf = CompactLFVals;

% Copy dimensions.
dim = [NumberSensors; NumberSourceSpaceNodes; NumberDirections];

% Close file.
fclose(FileID);
end
