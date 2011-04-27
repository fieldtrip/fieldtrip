function t = patientinfo(fid,seekset,t)
% T = PATIENTINFO(FID,SEEKSET,T);
%   FID = valid filepointer to a sqdfile
%   SEEKSET = starting point of file read
%   T = SQDHANDLE OBJECT/STRUCTURE
% Gets the info regarding patient (name, birthday,...) from the file pointed
% to by fid and returns a SQDHANDLE OBJECT if second argument is 
% passed as sqdhandle or if no arguemnt is passed or a structure is 
% passsed, it returns a structure


if nargin<1
    error('First argument must be valid file-pointer to a sqd-file');
elseif nargin<2
    seekset = -1;
end;

%********** Patient Information **********
% Get offset of patient information
fseek( fid, 32, seekset );
% long int patient_offset,patient_size,patient_maxcount,patient_count;
offset  = fread(fid,1,'long');
size    = fread(fid,1,'long');
maxcount= fread(fid,1,'long');
count   = fread(fid,1,'long');

% Read patient info
curoffset = offset;
curcount = 1;
tmp(curcount).id = '';
tmp(curcount).name = '';
tmp(curcount).birthdate = '';        
tmp(curcount).handedness = '';
while curoffset<offset+size*count
    fseek( fid, curoffset, seekset );
    infosize = fread(fid,1,'long');
    code = fread(fid,1,'int');
    subcode = fread(fid,1,'int');
    data = fread(fid,infosize,'uchar');
    dataend = min(find(data==0));
    data = setstr(data(1:dataend)');
    if subcode==1
        switch code
        case 1
            tmp(curcount).id = data;
        case 2
            tmp(curcount).name = data;
        case 3
            tmp(curcount).birthdate = data;
        case 4
            tmp(curcount).gender = deblank(data);
        case 5
            tmp(curcount).handedness = deblank(data);
        end;
    end;
    curoffset = curoffset+infosize;
    curcount = ceil((curoffset-offset)/size);
end;

if nargin>=2
    taglist = fieldnames(tmp);    
    if isa(t,'sqdhandle')
            t = set(t,'PatientInfo',tmp);
    elseif isstruct(t)
        for i =1:length(taglist)
            t = setfield(t,taglist{i},getfield(tmp,taglist{i}));
        end;
    end;
else
    t = tmp;
end;