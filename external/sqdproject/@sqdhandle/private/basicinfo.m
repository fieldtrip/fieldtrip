function t = basicinfo(fid,seekset,t)
% T = BASICINFO(FID,SEEKSET,T);
%   FID = valid filepointer to a sqdfile
%   SEEKSET = starting point of file read
%   T = SQDHANDLE OBJECT/STRUCTURE
% Gets the basic info from the file pointed to by fid and returns
% a SQDHANDLE OBJECT if second argument is passed as sqdhandle or
% if no arguemnt is passed or a structure is passsed, it returns a
% structure
if nargin<1
    error('First argument must be valid file-pointer to a sqd-file');
elseif nargin<2
    seekset = -1;
end;

% Get offset of basic information
fseek( fid, 16, seekset );
basic_offset = fread(fid,1,'long');
newline = sprintf('\n');
% Read basic information
fseek( fid, basic_offset, seekset );
tmp.Version	= fread(fid,1,'int');
tmp.Revision	= fread(fid,1,'int');
tmp.SystemID	= fread(fid,1,'int');
tmp.SystemName	= fread(fid,128,'uchar');
tmp.SystemName	= setstr(tmp.SystemName(1:...
    min(find(tmp.SystemName==0)))');
if ~isempty(tmp.SystemName)
    newindx = findstr(tmp.SystemName,newline);
    tmp.SystemName(newindx(end):end) = '';
end;
tmp.ModelName	= fread(fid,128,'uchar');
tmp.ModelName	= setstr(tmp.ModelName(1:...
    min(find(tmp.ModelName==0)))');
if ~isempty(tmp.ModelName)
    newindx = findstr(tmp.ModelName,newline);
    if size(newindx) > 0, tmp.ModelName(newindx(end):end) = ''; end;
end;
tmp.ChannelCount	= fread(fid,1,'int');
tmp.Comment	= fread(fid,256,'uchar');
tmp.Comment	= setstr(tmp.Comment'); 
if nargin>=2
    taglist = fieldnames(tmp);    
    if isa(t,'sqdhandle')
        for i =1:length(taglist)
            t = set(t,taglist{i},getfield(tmp,taglist{i}));
        end;
    elseif istruct(t)
        for i =1:length(taglist)
            t = setfield(t,taglist{i},getfield(tmp,taglist{i}));
        end;
    end;
else
    t = tmp;
end;