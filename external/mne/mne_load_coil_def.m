function [CoilDef,Header] = mne_load_coil_def(fname);
%
%
%   [CoilDef,Header] = mne_load_coil_def(fname);
%   CoilDef          = mne_load_coil_def(fname);
%
%   If file name is not specified, the standard coil definition file
%   $MNE_ROOT/setup/mne/coil_def.dat or $MNE_ROOT/share/mne/coil_def.dat is read
%
%   The content of the coil definition file is described in
%   section 5.6 of the MNE manual
%
%   This routine is modified from the original BrainStorm routine
%   created by John C. Mosher
%

%
%   John C. Mosher
%   Los Alamos National Laboratory
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Copyright (c) 2005 BrainStorm MMIV by the University of Southern California
%   Principal Investigator:
%   ** Professor Richard M. Leahy, USC Signal & Image Processing Institute
%
%
%   Revision 1.5  2006/09/08 19:27:13  msh
%   Added KIT coil type to mne_load_coil_def
%   Allow reading of measurement info by specifying just a file name.
%
%   Revision 1.4  2006/04/25 12:29:10  msh
%   Added Magnes and CTF reference coil drawings.
%
%   Revision 1.3  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.2  2006/04/17 15:01:34  msh
%   More small improvements.
%
%   Revision 1.1  2006/04/17 11:52:15  msh
%   Added coil definition stuff
%
%

me='MNE:mne_load_coil_def';

%% Setup the default inputs

if nargin == 0
    fname = mne_file_name('setup/mne/coil_def.dat');
    if ~exist(fname,'file')
        fname = mne_file_name('share/mne/coil_def.dat');
        if ~exist(fname,'file')
            error(me,'The standard coil definition file was not found');
        end
    end
end

% Read in the coil_def information

%% Open, Read in entire file, close

fid = fopen(fname,'rt');

if fid < 0
    error(me,'Could not open coil definition file %s',fname);
end

istr = 1; % string indexer
str_array = cell(1000,1); % preallocate
str_array{istr} = fgetl(fid); % get first string

while ischar(str_array{istr}),
    istr = istr + 1;
    str_array{istr} = fgetl(fid); % get next string
end
fclose(fid);

str_array = str_array(1:(istr-1)); % trim allocation

%% Read the Header, find the structure

% Gather the header lines
HeaderLines = strmatch('#',str_array); % lines that begin with a comment

% where are the structure lines
StructureLines = strmatch('# struct',str_array); % subset of header lines
% should be two
if length(StructureLines) ~= 2,
    error(me,'%s anticipates two lines of structures',mfilename)
end

% with each structure line is a format line
FormatLines = strmatch('# format',str_array); % subset of header lines
% assume there are also two

FieldNames = cell(1,2);
Format = cell(1,2);
% first structure is the coil information
% won't actually use the second structure, just its format
for i = 1:2,
    FieldNames{i}  = strread(str_array{StructureLines(i)},'%s');
    FieldNames{i}(1:2) = []; % strip the comment symbol and struct keyword
    
    [ignore,Format{i}] = strtok(str_array{FormatLines(i)},'''');
    Format{i} = strrep(Format{i},'''',''); % strip the single quotes
end

%% Allocate the arrays for loading

% interleave every fieldname with a null value
struct_arg = [FieldNames{1} cell(length(FieldNames{1}),1)]';
% each column an argument pair

% Preallocate a structure
[CoilDef(1:100)] = deal(struct(struct_arg{:}));


% Convert the rest of the string array to a structure
iCoil = 0; % counter
iLine = HeaderLines(end); % next line after header

while iLine < length(str_array), % number of lines in file
    iCoil = iCoil + 1; % next coil definition
    iLine = iLine + 1; % next line
    
    % first read the integer information on the coil
    % begin by breaking the line into two parts, numeric and description
    [numeric_items, description] = strtok(str_array{iLine},'"');
    temp = sscanf(numeric_items,Format{1}); % extra %s doesn't matter
    % assign temp in the order of the fields
    for i = 1:(length(FieldNames{1})-1),
        CoilDef(iCoil).(FieldNames{1}{i}) = temp(i);
    end
    % then assign the description
    % let's strip the quotes first
    description = strrep(description,'"','');
    CoilDef(iCoil).(FieldNames{1}{end}) = description;
    
    % now read the coil definition
    CoilDef(iCoil).coildefs = zeros(CoilDef(iCoil).num_points,7);
    for i = 1:CoilDef(iCoil).num_points,
        iLine = iLine + 1;
        CoilDef(iCoil).coildefs(i,:) = sscanf(str_array{iLine},...
            Format{2})';
    end
    
    % now draw it
    % local subfunction below
    CoilDef(iCoil).FV = draw_sensor(CoilDef(iCoil));
    
end

CoilDef = CoilDef(1:iCoil); % trim allocation

Header = str_array(HeaderLines);

function FV = draw_sensor(CoilDef);

% create a patch based on the sensor type

% The definitions as of 14 October 2005:
% for i = 1:3:length(CoilDef),fprintf('%d %s\n',CoilDef(i).id,CoilDef(i).description);end
% 2 Neuromag-122 planar gradiometer size = 27.89  mm base = 16.20  mm
% 2000 Point magnetometer
% 3012 Vectorview planar gradiometer T1 size = 26.39  mm base = 16.80  mm
% 3013 Vectorview planar gradiometer T2 size = 26.39  mm base = 16.80  mm
% 3022 Vectorview magnetometer T1 size = 25.80  mm
% 3023 Vectorview magnetometer T2 size = 25.80  mm
% 3024 Vectorview magnetometer T3 size = 21.00  mm
% 4001 Magnes WH2500 magnetometer size = 11.50  mm
% 4002 Magnes WH3600 gradiometer size = 18.00  mm base = 50.00  mm
% 5001 CTF axial gradiometer size = 18.00  mm base = 50.00  mm
% 7001 BabySQUID I axial gradiometer size = 6.0 mm base = 50 mm
% 7002 BabyMEG magnetometer 
% 7003 BabyMEG reference magnetometer 

FV = struct('faces',[],'vertices',[]); % standard convention

% recall that vertices are 1 per ROW, not column, of matrix

switch CoilDef.id
    case {2,3012,3013,3011}
        % square figure eight
        % wound by right hand rule such that +x side is "up" (+z)
        LongSide = CoilDef.size*1000; % length of long side in mm
        Offset = 2.5; % mm offset of the center portion of planar grad coil
        FV.vertices = [0 0 0; Offset 0 0; ...
            Offset -LongSide/2 0; LongSide/2 -LongSide/2 0; ...
            LongSide/2 LongSide/2 0; ...
            Offset LongSide/2 0; Offset 0 0; ...
            0 0 0; -Offset 0 0; -Offset -LongSide/2 0; ...
            -LongSide/2 -LongSide/2 0; ...
            -LongSide/2 LongSide/2 0; ...
            -Offset LongSide/2 0; -Offset 0 0]/1000;
        FV.faces = [1:length(FV.vertices)];
    case 2000
        % point source
        LongSide = 2; % mm, tiny square
        FV.vertices = [-1 1 0;1 1 0;1 -1 0; -1 -1 0]*LongSide/1000/2;
        FV.faces = [1:length(FV.vertices)];
    case {3022, 3023, 3024}
        % square magnetometer
        LongSide = CoilDef.size*1000; % mm, length of one side
        FV.vertices = [-1 1 0;1 1 0;1 -1 0; -1 -1 0]*LongSide/1000/2;
        FV.faces = [1:length(FV.vertices)];
    case {4001,4003,5002,7002,7003}
        % round magnetometer
        Radius = CoilDef.size*1000/2; % mm, radius of coil
        Len_cir = 15; % number of points for circle
        circle = cos(2*pi*[0:(Len_cir-1)]/Len_cir) + ...
            sqrt(-1)*sin(2*pi*[0:(Len_cir-1)]/Len_cir); % complex circle unit
        FV.vertices = ...
            [real(circle)' imag(circle)' zeros(Len_cir,1)]*Radius/1000;
        FV.faces = [1:length(FV.vertices)];
    case {4002, 5001, 5003, 4004, 6001, 7001}
        % round coil 1st order gradiometer
        Radius = CoilDef.size*1000/2; % mm radius
        Baseline = CoilDef.baseline*1000; % axial separation
        Len_cir = 15; % number of points for circle
        % This time, go all the way around circle to close it fully
        circle = cos(2*pi*[0:Len_cir]/Len_cir) + ...
            sqrt(-1)*sin(2*pi*[0:Len_cir]/Len_cir); % complex circle unit
        circle = circle*Radius; % scaled
        FV.vertices = ...
            [[real(circle)' imag(circle)' zeros(Len_cir+1,1)];... % first coil
            [real(circle)' -imag(circle)' zeros(Len_cir+1,1)+Baseline]]/1000; % 2nd coil
        FV.faces = [1:length(FV.vertices)];
    case {5004,4005}
        % round coil 1st order off-diagonal gradiometer
        Radius = CoilDef.size*1000/2; % mm radius
        Baseline = CoilDef.baseline*1000; % axial separation
        Len_cir = 15; % number of points for circle
        % This time, go all the way around circle to close it fully
        circle = cos(2*pi*[0:Len_cir]/Len_cir) + ...
            sqrt(-1)*sin(2*pi*[0:Len_cir]/Len_cir); % complex circle unit
        circle = circle*Radius; % scaled
        FV.vertices = ...
            [[real(circle)'+Baseline/2.0 imag(circle)' zeros(Len_cir+1,1)];... % first coil
            [real(circle)'-Baseline/2.0 -imag(circle)' zeros(Len_cir+1,1)]]/1000; % 2nd coil
        FV.faces = [1:length(FV.vertices)];
    otherwise
        FV = [];
end


