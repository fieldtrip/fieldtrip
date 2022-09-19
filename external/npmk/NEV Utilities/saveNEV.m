function saveNEV(NEV, varargin)

%% 
% Save an .NEV file from an NEV structure (gained by using openNEV)
% Works with file spec 2.3
% 
% Use saveNEV(NEV, filename, noreport)

% All input arguments are optional. Input arguments can be in any order.
%
%   NEV:        Name of the NEV structure to be saved.
%               DEFAULT: The workspace NEV structure in the workspace will
%               be saved.
%
%   filename:   Name of the new NEV file name.
%               DEFAULT: -out.nev will be added to the end of the current
%               file name.
%
%   'noreport': Will not display status reports or warnings.
%               DEFAULT: will display status reports and warnings.
%
%
%   OUTPUT:     None
%
%   Example 1: 
%   saveNEV
%
%   In the example above, a new file containing the NEV structure in the
%   workspace will be saved. The file name will have -out.NEV added to its
%   name. Statuses of the progress of saveNEV and also warnings about the
%   risks of saving will be displayed.
%
%   Example 2:
%   openNSx(NEV, 'myNewNEVFile.nev', 'noreport');
%
%   In the example above, a new file containing the NEV structure in the
%   workspace will be saved. The file name will be myNewNEVFile.nev. 
%   Statuses of the progress of saveNEV and also warnings about the
%   risks of saving will be not be displayed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Nick Halper
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version History
%
% 1.0.0.0: Nick Halper
%   - Initial release.
%
% 1.1.0.0: Kian Torab
%   - Added the ability to suppress the fiel saving warning with
%     'noreport' input parameter.
%   - Re-structured help to better match the suite's style. Added examples.
%   - Added the 'noreport' key to supress statuses and warnings.
%   - Improved error checking for input arguments.
%
% 1.2.0.0: Nick Halper - 19/8/29
%   - Allows saveNEV to dynamically determine header offset and total
%   number of extended headers instead of relying on info in MetaTags/Raw
%   Data of NEV structure.
%   - Fixed a bug where choosing to overwrite would crash the script
%   - Fixed an issue where bank letters were not being written correctly, writing
%   blank values
%   - Modified the script to work with 256 channel files
%
% 1.3.0.0: Stephen hou - 19/8/30
%   - Implemented saving of files that contain NeuroMotive/tracking data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verify FilePath and establish overwrite paramaters
if not(NEV.MetaTags.FileSpec == '2.3')
    disp(strcat('This function only functions on file spec 2.3,;this is your file spec:',NEV.MetaTags.FileSpec));
    return
end

if nargin > 1
    for idx = 1:nargin-1
        if strcmpi(varargin{idx}, 'noreport')
            reportFlag = 0;
        elseif length(varargin{idx})>3 && ...
                (strcmpi(varargin{idx}(3),'\') || ...
                 strcmpi(varargin{idx}(1),'/') || ...
                 strcmpi(varargin{idx}(2),'/') || ...
                 strcmpi(varargin{idx}(1:2), '\\')) 
            FilePath = varargin{1};
        end
    end
end

% Validating input arguments
if ~exist('FilePath',   'var'); FilePath = [fullfile(NEV.MetaTags.FilePath,NEV.MetaTags.Filename) '-out.nev']; end
if ~exist('reportFlag', 'var'); reportFlag = 1; end;

% Warning user about the consequences of the modified NEV file
if reportFlag
    Accept = input('This script will save a new file with the .NEV extensions, but you should retain the previous file. Do you acknowledge the risk inherent in saving modified versions of data files? (Y/N)','s');
    if ~strcmpi(Accept,'y')
        disp('Ending Script...');
        return
    end
end

% Validating and opening the writable file
if exist(FilePath)
    
        if exist(FilePath)
         
                disp('File already exists!');
                OverwritePrompt = input('Would you like to overwrite? (Y/N)','s');
                if strcmpi(OverwritePrompt,'y')
                    Overwrite = 1;
                    delete(FilePath);
                else
                    while exist(FilePath)
                        ExistCount = ExistCount + 1;
                        FilePath = fullfile(NEV.MetaTags.FilePath,[NEV.MetaTags.Filename '-Aligned-',num2str(ExistCount),NEV.MetaTags.FileExt]);
                    end
                    
                end
            
        end
end

FileID = fopen(FilePath, 'w', 'ieee-le');
    
if (FileID <= 0)
    disp('No file was selected.');
    return;
end

%% Write the basic header into the file

% General warining about the length of time it takes to save a NEV file so
% the user does not abort the process
disp(['Saving NEV file ' FilePath '. This can take a long time. Please be patient.']);

if reportFlag; disp('Writing Basic Header...'); end

fwrite(FileID,NEV.MetaTags.FileTypeID(1:8));
fwrite(FileID, [str2double(NEV.MetaTags.FileSpec(1)) str2double(NEV.MetaTags.FileSpec(3))], 'uint8');
fwrite(FileID,str2double(NEV.MetaTags.Flags),'uint16');
fwrite(FileID,NEV.MetaTags.HeaderOffset,'uint32');
fwrite(FileID,NEV.MetaTags.PacketBytes,'uint32');
fwrite(FileID,NEV.MetaTags.SampleRes,'uint32');
fwrite(FileID,NEV.MetaTags.TimeRes,'uint32');
fwrite(FileID,NEV.MetaTags.DateTimeRaw,'uint16');
fwrite(FileID,'saveNEV$version1001$$$$$$$$$$$$$'); 
fwrite(FileID,NEV.MetaTags.Comment); 
ExtendedHeaderBytes = NEV.MetaTags.HeaderOffset-ftell(FileID)+4;
fwrite(FileID,ExtendedHeaderBytes/32,'uint32');

EndOfBasicHeader = ftell(FileID);

%%
% Write the extended header into the file. 

%Handling packets with array information
if isfield(NEV,'ArrayInfo')
    if reportFlag; disp('Writing Array Header...'); end
    if isfield(NEV.ArrayInfo,'ElectrodeName')
        fwrite(FileID,'ARRAYNME');
        fwrite(FileID,NEV.ArrayInfo.ElectrodeName); %Must null terminate
    end
    if isfield(NEV.ArrayInfo,'ArrayComment')
        fwrite(FileID,'ECOMMENT');
        fwrite(FileID,NEV.ArrayInfo.ArrayComment); %Must null terminate
    end
    if isfield(NEV.ArrayInfo,'ArrayCommentCont')
        fwrite(FileID,'CCOMMENT');
        fwrite(FileID,NEV.ArrayInfo.ArrayCommentCont); %Must null terminate
    end
    if isfield(NEV.ArrayInfo,'MapFile')
        fwrite(FileID,'MAPFILE'); %+NULL
        fwrite(FileID,NEV.ArrayInfo.MapFile); %Must null terminate
    end
end

if isfield(NEV,'ElectrodesInfo')
    if reportFlag; disp('Writing Electrode Header...'); end
    if (isfield(NEV.ElectrodesInfo(1),'ElectrodeID'))
        
    %Find length of electrode count, loop through for that count and fill
    %in  NEUEVWAV packets. 
        for IDX = 1:length(NEV.ElectrodesInfo)
            Before = ftell(FileID);
            fwrite(FileID,'NEUEVWAV');
            fwrite(FileID,NEV.ElectrodesInfo(IDX).ElectrodeID,'uint16');
            

            %fwrite(FileID,NEV.ElectrodesInfo(IDX).ConnectorBank);
            switch NEV.ElectrodesInfo(IDX).ConnectorBank
                case 'A'
                    fwrite(FileID,1,'uint8');
                case 'B'
                    fwrite(FileID,2,'uint8');
                case 'C'
                    fwrite(FileID,3,'uint8');
                case 'D'
                    fwrite(FileID,4,'uint8');
                case 'E'
                    fwrite(FileID,5,'uint8');
                case 'F'
                    fwrite(FileID,6,'uint8');
                case 'G'
                    fwrite(FileID,7,'uint8');
                case 'H'
                    fwrite(FileID,8,'uint8');
                case 'I'
                    fwrite(FileID,9,'uint8');
            end
        
            fwrite(FileID,NEV.ElectrodesInfo(IDX).ConnectorPin,'uint8');
            fwrite(FileID,NEV.ElectrodesInfo(IDX).DigitalFactor,'uint16');
            fwrite(FileID,NEV.ElectrodesInfo(IDX).EnergyThreshold,'uint16');
            fwrite(FileID,NEV.ElectrodesInfo(IDX).HighThreshold,'int16');
            fwrite(FileID,NEV.ElectrodesInfo(IDX).LowThreshold,'int16');
            fwrite(FileID,NEV.ElectrodesInfo(IDX).Units,'uint8');
            fwrite(FileID,NEV.ElectrodesInfo(IDX).WaveformBytes,'uint8');
            if isempty(NEV.Data.Spikes.Waveform)
                fwrite(FileID,48,'uint16');
                SpikeLength = 48;
            else
                fwrite(FileID,length(NEV.Data.Spikes.Waveform(:,1)),'uint16');
                SpikeLength = length(NEV.Data.Spikes.Waveform(:,1));
            end
            %if file type is 2.2, don't need previous field and end in 10
            %zeros
            fwrite(FileID,zeros(8,1),'uint8');
            After = ftell(FileID);
            if After-Before ~= 32
                disp('Broken Electrode Info')
                NEV.ElectrodesInfo(IDX).ConnectorBank
            end
        end
    end
    if (isfield(NEV.ElectrodesInfo(1),'ElectrodeLabel'))
        for IDX = 1:length(NEV.ElectrodesInfo)
        Before = ftell(FileID);
        fwrite(FileID,'NEUEVLBL');
        fwrite(FileID, NEV.ElectrodesInfo(IDX).ElectrodeID,'uint16');
        fwrite(FileID, NEV.ElectrodesInfo(IDX).ElectrodeLabel);%Must be nulll terminated
        fwrite(FileID, zeros(6,1),'uint8');
        After = ftell(FileID);
        if After-Before ~= 32
            disp('Broken Electrode Label')
        end
        end
    end
    if (isfield(NEV.ElectrodesInfo(1),'HighFreqCorner'))
        for IDX = 1:length(NEV.ElectrodesInfo)
        Before = ftell(FileID);
        fwrite(FileID,'NEUEVFLT');
        fwrite(FileID, NEV.ElectrodesInfo(IDX).ElectrodeID,'uint16');
        fwrite(FileID, NEV.ElectrodesInfo(IDX).HighFreqCorner,'uint32');
        fwrite(FileID, NEV.ElectrodesInfo(IDX).HighFreqOrder,'uint32');
        fwrite(FileID, NEV.ElectrodesInfo(IDX).HighFilterType,'uint16');
        fwrite(FileID, NEV.ElectrodesInfo(IDX).LowFreqCorner,'uint32');
        fwrite(FileID, NEV.ElectrodesInfo(IDX).LowFreqOrder,'uint32');
        fwrite(FileID, NEV.ElectrodesInfo(IDX).LowFilterType,'uint16');
        fwrite(FileID, zeros(2,1), 'uint8');
        After = ftell(FileID);
        if After-Before ~= 32
            disp('Broken High Freq Corner');
        end
        end
    end   
end
%Digital inputs
if isfield(NEV,'IOLabels')
    if reportFlag; disp('Writing IOLabels Header...'); end
        for IDX = 1:length(NEV.IOLabels)
            Before = ftell(FileID);
            fwrite(FileID,'DIGLABEL');
            fwrite(FileID,'Serial1XXXXXXXXX','uint8');
            fwrite(FileID, IDX - 1, 'uint8');
            fwrite(FileID, zeros(7,1),'uint8');
            After = ftell(FileID);
            if After-Before ~= 32
               disp('Broken IO Labels');
            end
        end
end
    
%Video Packets
if isfield(NEV,'VideoSyncInfo')
    if reportFlag; disp('Writing Video Header...'); end
    for IDX = 1:length(NEV.VideoSyncInfo)
    Before = ftell(FileID);
    fwrite(FileID,'VIDEOSYN');
    fwrite(FileID, NEV.VideoSyncInfo(IDX).SourceID, 'uint16');
    fwrite(FileID, NEV.VideoSyncInfo(IDX).SourceName(1:16));
    fwrite(FileID, NEV.VideoSyncInfo(IDX).FrameRateFPS,'single');
    fwrite(FileID, zeros(2,1),'uint8');
    After = ftell(FileID);
    if After-Before ~= 32
       disp('Broken Video Sync Info');
       PacketNumber = IDX;
       TotalPackets = length(NEV.VideoSyncInfo);
    end
    end
end

if isfield(NEV,'NSAS')
    %This might exist in a future version of Central
end

if isfield(NEV,'ObjTrackInfo')
    if reportFlag; disp('Writing Tracking Header...'); end
    for IDX = 1:length(NEV.ObjTrackInfo)
    Before = ftell(FileID);
    fwrite(FileID,'TRACKOBJ');
    fwrite(FileID, NEV.ObjTrackInfo(IDX).TrackableType,'uint16');
    fwrite(FileID, NEV.ObjTrackInfo(IDX).TrackableID,'uint32');%This is an error and should be two different uint16 values, but we can read it back into file this way.
    NEV.ObjTrackInfo(IDX).TrackableName = pad(NEV.ObjTrackInfo(IDX).TrackableName,16,'X'); 
    fwrite(FileID, NEV.ObjTrackInfo(IDX).TrackableName);
    fwrite(FileID, zeros(2,1),'uint8');
    After = ftell(FileID);
    if After-Before ~= 32
       disp('Broken Obj Track Info');
    end
    end
end

if isfield(NEV,'Rabbits')
    %Fill in the details about Rabbits at some point in the future.
end

EndOfExtendedHeader = ftell(FileID);

%% Edit the Basic Header to Account for Number of Extended Headers

fseek(FileID,12,'bof');
fwrite(FileID,EndOfExtendedHeader,'uint32');
fseek(FileID,322,'bof');
fwrite(FileID,(EndOfExtendedHeader-EndOfBasicHeader)/32,'uint32');
fseek(FileID,EndOfExtendedHeader,'bof');


%%
% Write Data
BytesInPackets = NEV.MetaTags.PacketBytes;
Broken = 0;

%SerialDigitalIO CHECK
if ~isempty(NEV.Data.SerialDigitalIO.TimeStamp)
    if reportFlag; disp('Writing Serial/Digital Data...'); end
    for IDX = 1:length(NEV.Data.SerialDigitalIO.TimeStamp)
        Before = ftell(FileID);
        fwrite(FileID, NEV.Data.SerialDigitalIO.TimeStamp(IDX),'uint32');
        %ftell(FileID)-Before
        fwrite(FileID, 0,'uint16');
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.SerialDigitalIO.InsertionReason(IDX));
        %ftell(FileID)-Before
        fwrite(FileID, '0');
        %ftell(FileID)-Before
        if ~isempty(NEV.Data.SerialDigitalIO.Value)
            fwrite(FileID, NEV.Data.SerialDigitalIO.Value(IDX),'uint16');
        else
            fwrite(FileID, NEV.Data.SerialDigitalIO.UnparsedData(IDX),'uint16');
        end
        %ftell(FileID)-Before
        fwrite(FileID, zeros(BytesInPackets-10,1),'uint8');
        %ftell(FileID)-Before
        After = ftell(FileID);
        if After-Before ~= BytesInPackets
            Broken = 1;
            %After-Before
            %CurrentPacket = IDX
            %TotalPackets = length(NEV.Data.SerialDigitalIO.TimeStamp)
        end
    end
    if Broken == 1
        disp('Serial Digital Packet Corrupted');
        Broken = 0;
    end
end

%Spikes CHECK
if ~isempty(NEV.Data.Spikes.TimeStamp)
    if reportFlag; disp('Writing Spike Data...'); end
    for IDX = 1:length(NEV.Data.Spikes.TimeStamp)
        Before = ftell(FileID);
        fwrite(FileID, NEV.Data.Spikes.TimeStamp(IDX),'uint32');

        fwrite(FileID, NEV.Data.Spikes.Electrode(IDX),'uint16');

        fwrite(FileID, NEV.Data.Spikes.Unit(IDX),'uchar');

        fwrite(FileID, 0,'uchar');

        fwrite(FileID, NEV.Data.Spikes.Waveform(:,IDX)','int16');
        
        
        
        %for Value = 1:SpikeLength
        %    fwrite(FileID, NEV.Data.Spikes.Waveform(Value,IDX),'int16');
        %end
        
        After = ftell(FileID);
        if After-Before ~= BytesInPackets
            Broken = 1;
        end
    end
    if Broken == 1
        disp('Spike Packet Corrupted')
        Broken = 0;
    end
end
%disp('done')
%Comments CHECK
if ~isempty(NEV.Data.Comments.TimeStamp)
    if reportFlag; disp('Writing Comment Data...'); end
    for IDX = 1:length(NEV.Data.Comments.TimeStamp)
        Before = ftell(FileID);
        fwrite(FileID, NEV.Data.Comments.TimeStamp(IDX),'uint32');
        %ftell(FileID)-Before
        fwrite(FileID, 65535, 'uint16');
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.Comments.CharSet(IDX),'uint8');
        %ftell(FileID)-Before
        fwrite(FileID, 0,'uint8');
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.Comments.Color(IDX),'uint32');
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.Comments.Text(IDX,:));
        %ftell(FileID)-Before
        %Need to handle extra characters here etc
        fwrite(FileID, zeros(BytesInPackets-(ftell(FileID)-Before),1),'uint8');
        After = ftell(FileID);
        if After-Before ~= BytesInPackets
            Broken = 1;
            %After-Before
            %CurrentPacket = IDX
            %TotalPackets = length(NEV.Data.Comments.TimeStamp)
        end
    end
    if Broken == 1
        disp('Comment Packet Corrupted')
        Broken = 0;
    end
end

if ~isempty(NEV.Data.TrackingEvents.TimeStamp)
    if reportFlag; disp('Writing Tracking Event Data'); end
    for IDX = 1:length(NEV.Data.TrackingEvents.TimeStamp)
        Before = ftell(FileID);
        fwrite(FileID, NEV.Data.TrackingEvents.TimeStamp(IDX),'uint32');
        fwrite(FileID, 65535,'uint16');
        fwrite(FileID, 255, 'uint8');
        %fwrite(FileID,  2, 'uint8');
        fwrite(FileID, 0, 'uint8'); %placeholder until I ask hyrum what the Nreuromotive flag is
        fwrite(FileID, int2str(NEV.Data.TrackingEvents.ROINum(IDX)),'uint8');
        %byte 1 is ROI #, byte 2 is enter or exit -- is it 1 and 2?
        if strcmp(NEV.Data.TrackingEvents.Event(IDX),'Enter')
            fwrite(FileID, 1,'uint8');
        elseif strcmp(NEV.Data.TrackingEvents.Event(IDX),'Exit')
            fwrite(FileID, 2,'uint8');
        end
        fwrite(FileID, zeros(2,1), 'uint8');
        fwrite(FileID, strcat(NEV.Data.TrackingEvents.ROIName{IDX},':',int2str(NEV.Data.TrackingEvents.ROINum(IDX)),':',NEV.Data.TrackingEvents.Event{IDX},':',int2str(NEV.Data.TrackingEvents.Frame(IDX)),':','placeholder until i figure out wtf this is'));
        fwrite(FileID, zeros(BytesInPackets-(ftell(FileID)-Before),1),'uint8');
        After = ftell(FileID);
        if After-Before ~= BytesInPackets
            Broken = 1;
            %After-Before
            %CurrentPacket = IDX
            %TotalPackets = length(NEV.Data.Comments.TimeStamp)
        end
    end
    if Broken == 1
        disp('Tracking Event Packet Corrupted')
        Broken = 0;
    end
end

if ~isempty(NEV.Data.VideoSync.TimeStamp)
    if reportFlag; disp('Writing VideoSync Data...'); end
    for IDX = 1:length(NEV.Data.VideoSync.TimeStamp)
        Before = ftell(FileID);
        fwrite(FileID, NEV.Data.VideoSync.TimeStamp(IDX),'uint32');
        %ftell(FileID)-Before
        fwrite(FileID, 65534, 'uint16');
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.VideoSync.FileNumber(IDX),'uint16');
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.VideoSync.FrameNumber(IDX),'uint32');%Wrong Size
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.VideoSync.ElapsedTime(IDX),'uint32');%Wrong Size
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.VideoSync.SourceID(IDX),'uint32');
        %ftell(FileID)-Before
        fwrite(FileID, zeros(BytesInPackets - 20,1),'uint8');
        %ftell(FileID)-Before
        After = ftell(FileID);
        if After-Before ~= BytesInPackets
            Broken = 1;
            %After-Before
            %CurrentPacket = IDX
            %TotalPackets = length(NEV.Data.VideoSync.TimeStamp)
        end
        
    end
    if Broken == 1
        disp('Video Sync Packet Corrupted')
        Broken = 0;
    end
end

if ~isempty(NEV.Data.Tracking)
    if reportFlag; disp('Writing Tracking Data...'); end
    TrackingFieldNames = fieldnames(NEV.Data.Tracking);
        for TrackingField = 1:numel(TrackingFieldNames)
            for IDX = 1:length(NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).TimeStamp)
                Before = ftell(FileID);
                if NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).MarkerCount == 0
                    
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).TimeStamp(IDX),'uint32');
                    %ftell(FileID)-Before
                    fwrite(FileID, 65533, 'uint16');
                    %ftell(FileID)-Before
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).ParentID(IDX),'uint16');
                    %ftell(FileID)-Before
                    fwrite(FileID, TrackingField-1,'uint16'); %Node ID
                    %ftell(FileID)-Before                    
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).NodeCount(IDX),'uint16');
                    %ftell(FileID)-Before
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).MarkerCount(IDX),'uint16');
                    %ftell(FileID)-Before
                    fwrite(FileID, zeros(2,1),'uint8');
                    %ftell(FileID)-Before
                    fwrite(FileID, zeros(BytesInPackets - 16,1),'uint8');
      
                else
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).TimeStamp(IDX),'uint32');
                    %ftell(FileID)-Before
                    fwrite(FileID, 65533, 'uint16');
                    %ftell(FileID)-Before
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).ParentID(IDX),'uint16');
                    %ftell(FileID)-Before
                    fwrite(FileID, TrackingField-1,'uint16'); %Node ID
                    %ftell(FileID)-Before
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).NodeCount(IDX),'uint16');
                    %ftell(FileID)-Before
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).MarkerCount(IDX),'uint16');
                    %ftell(FileID)-Before
                    MarkerCoordinatesBytes = length(NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).MarkerCoordinates(IDX).X);
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).MarkerCoordinates(IDX).X,'uint16');
                    %ftell(FileID)-Before
                    fwrite(FileID, NEV.Data.Tracking.(TrackingFieldNames{TrackingField}).MarkerCoordinates(IDX).Y,'uint16');

                    %ftell(FileID)-Before
                    if MarkerCoordinatesBytes == 1
                        fwrite(FileID, zeros(BytesInPackets - 14-MarkerCoordinatesBytes*4,1),'uint8');
                    else
                        fwrite(FileID, zeros(BytesInPackets - 14-MarkerCoordinatesBytes*4,1),'uint8');
                    end
                end
                    
                After = ftell(FileID);
                if After-Before ~= BytesInPackets
                    Broken = 1;
                    %After-Before
                    %CurrentPacket = IDX
                    %TotalPackets = length(NEV.Data.VideoSync.TimeStamp)
                end
                %Must somehow terminate in correct number of zeros
            end
            if Broken == 1
                disp('Tracking Packet Corrupted')
                disp(After - Before)
                Broken = 0;
            end
        end
end

if ~isempty(NEV.Data.PatientTrigger.TimeStamp)
    if reportFlag; disp('Writing Patient Trigger Data...'); end
    for IDX = 1:length(NEV.Data.PatientTrigger.TimeStamp)
        Before = ftell(FileID);
        fwrite(FileID, NEV.Data.PatientTrigger.TimeStamp(IDX),'uint32');
        %ftell(FileID)-Before
        fwrite(FileID, 65532, 'uint16');
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.PatientTrigger.TriggerType(IDX),'uint16');
        %ftell(FileID)-Before
        fwrite(FileID, zeros(BytesInPackets - 8, 1),'uint8');
        %ftell(FileID)-Before
        After = ftell(FileID);
        if After-Before ~= BytesInPackets
            Broken = 1;
            %After-Before
            %CurrentPacket = IDX
            %TotalPackets = length(NEV.Data.PatientTrigger.TimeStamp)
        end
    end
    if Broken == 1
        disp('Patient Trigger Packet Corrupted')
        Broken = 0;
    end
    
end

if ~isempty(NEV.Data.Reconfig.TimeStamp)
    if reportFlag; disp('Writing Reconfig Data...'); end
    for IDX = 1:length(NEV.Data.Reconfig.TimeStamp)
        Before = ftell(FileID);
        fwrite(FileID, NEV.Data.Reconfig.TimeStamp(IDX),'uint32');
        %ftell(FileID)-Before
        fwrite(FileID, 65531, 'uint16');
        %ftell(FileID)-Before
        fwrite(FileID, NEV.Data.Reconfig.ChangeType(IDX),'uint16');
        %ftell(FileID)-Before
        fwrite(FileID, zeros(BytesInPackets - 8,1),'uint8');
        %ftell(FileID)-Before
        After = ftell(FileID);
        if After-Before ~= BytesInPackets
            Broken = 1;
            %After-Before
            %CurrentPacket = IDX
            %TotalPackets = length(NEV.Data.Reconfig.TimeStamp)
        end
    end
    if Broken == 1
        disp('Reconfig Packet Corrupted')
        Broken = 0;
    end
end

if reportFlag; disp('Finished!'); end

clear After
clear Before
clear Broken
clear BytesInPackeets
clear ExtendedHederBytes
clear FilePath
clear IDX
clear SpikeLength
clear Value

fclose('all');

