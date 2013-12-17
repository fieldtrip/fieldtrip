function [dim lf] = readBESAlft(filename)

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
% Return:
%     [dim] 
%         Leadfield dimensions {<number sensors>, <number source space
%         nodes>, <number source directions>}
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

% Read compactly stored LF values (int16)
NumberLFValues = NumberColumns * NumberSensors;
[CompactLFVals NrReadElements] = fread(FileID, NumberLFValues, 'int16');
if(NrReadElements ~= NumberLFValues)
	fprintf('Could not read leadfield values from file.\n');
	return
end

% Reshape matrix with compact LF values.
CompactLFVals = reshape(CompactLFVals, NumberColumns, NumberSensors);
CompactLFVals = CompactLFVals';

% Compute factors for converting compactly stored LF values
% to float values.
ConversionFactors = MaxVals / 32767.0;
clear MaxVals;

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
