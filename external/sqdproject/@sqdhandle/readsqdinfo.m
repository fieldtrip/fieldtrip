function t = readsqdinfo(t,fname)
% READSQDINFO internal function which reads the given sqd-file

seekset=-1;

if nargin<2|~ischar(fname),error('Incorrect arguments');end;
if isempty(fname)       % Default params
    t.Version       = 1;
    t.Revision      = 1;
    t.SystemID      = '';
    t.SystemName    = '';
    t.ModelName     = '';
    t.FileName      = '';
    t.Amplifier     = [];
    t.ChannelCount	= 192;
    t.Channels      = chanhandle(t.FileName,0);
    t.Comment       = '';
    t.PatientInfo   = [];
    acqtype	        = acqparam;
else    
    % open file as binary and little endian
    fid = fopen(fname,'rb','l');
    if fid==-1,
        disp('File not found');
        [f,p] = uigetfile('*.sqd','MEG150 files');
        if ((~isstr(f)) | ~min(size(f))), 
            t = [];
            return; 
        end;
        fname = [p f];
        t.FileName = [p f];
        fid = fopen(fname,'rb','l');
    end;
    % Read in basic info - verion, revision, channel count, etc.
    t = basicinfo(fid,seekset,t);
    % Read in patient info
    t = patientinfo(fid,seekset,t);
    % Read in Amplifier Gain values
    t.Amplifier     = amplifierinfo(fid,seekset);
    % Read acquisition parameters
    t.acqparam = acqparam(fid);
    % Close file
    fclose(fid);
end;
% Update channel info
for channum = 0:t.ChannelCount-1
    t.Channels(channum+1) = chanhandle(t.FileName,channum);
end;
