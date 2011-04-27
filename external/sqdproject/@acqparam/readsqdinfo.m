function t = readsqdinfo(t,fileinfo,acq_offset,raw_offset)
% t = readsqdinfo(t,fileinfo,acq_offset,raw_offset); Read info from sqd-file
% and write it to ACQPARAM object

SEEK_SET=-1;    % Start of the file for fseek: -1 => at the beginning

if ischar(fileinfo) % If filename is passed
    % open file as binary and little endian
    fid = fopen(fileinfo,'rb','l');
    if fid==-1,
        disp('File not found');
        [f,p] = uigetfile('*.sqd','MEG160 files');
        if ((~isstr(f)) | ~min(size(f))), return; end;
        fname = [p f];
        fid = fopen(fname,'rb','l');
    end;
elseif isnumeric(fileinfo)% If filepointer is passed
    if fileinfo>0
        fid = fileinfo;
    else
        error('Invalid fid');
    end;
end;

% Initialize the offsets
switch nargin
case 2
    % Get offset of acquisition parameters
    fseek( fid, 128, SEEK_SET );
    acq_offset = fread(fid,1,'long');
    % Get Acquisition type
    fseek(fid,acq_offset,SEEK_SET);
    t.AcquisitionType       = fread(fid,1,'long');
    % Get Raw Offset
    fseek(fid,t.Private.Offset2Raw(t.AcquisitionType),SEEK_SET);
    raw_offset  = fread(fid,1,'long');
case 3
    % Get Acquisition type
    fseek(fid,acq_offset,SEEK_SET);
    t.AcquisitionType       = fread(fid,1,'long');
    % Get Raw Offset
    fseek(fid,t.Private.Offset2Raw(t.AcquisitionType),SEEK_SET);
    raw_offset  = fread(fid,1,'long');
end;
%Save Offsets
t.RawOffset         = raw_offset;

% Read acquisition type
fseek( fid, acq_offset, SEEK_SET );
t.AcquisitionType       = fread(fid,1,'long');
t.SampleRate            = fread(fid,1,'double');
t.Datatype              = t.Private.Datatypes{t.AcquisitionType};
 
switch t.AcquisitionType
    case 1  
        % Read acquisition parameters
        t.SampleInfo.SamplesAcquired    = fread(fid,1,'long');
        t.SampleInfo.ActSamplesAcquired = fread(fid,1,'long');
        t.SampleInfo.SamplesAvailable   = t.SampleInfo.ActSamplesAcquired;
    case {2,3}
        % Read acquisition parameters
        t.SampleInfo.FrameLength        = fread(fid,1,'long');
        t.SampleInfo.SamplesAvailable   = t.SampleInfo.FrameLength;
        t.SampleInfo.PreTrigger         = fread(fid,1,'long');
        t.SampleInfo.AverageCount       = fread(fid,1,'long');
        t.SampleInfo.ActAverageCount    = fread(fid,1,'long');
end;
if ischar(fileinfo);fclose(fid);end;
return;
