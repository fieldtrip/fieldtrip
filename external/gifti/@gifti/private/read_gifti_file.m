function this = read_gifti_file(filename, this)
% Low level reader of GIfTI 1.0 files
% FORMAT this = read_gifti_file(filename, this)
% filename    - XML GIfTI filename
% this        - structure with fields 'metaData', 'label' and 'data'.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id$

% Import XML-based GIfTI file
%--------------------------------------------------------------------------
try
    t = xmltree(filename);
catch
    error('[GIFTI] Loading of XML file %s failed.', filename);
end

% Root element of a GIFTI file
%--------------------------------------------------------------------------
if ~strcmp(get(t,root(t),'name'),'GIFTI')
    error('[GIFTI] %s is not a GIFTI 1.0 file.', filename);
end
attr = cell2mat(attributes(t,'get',root(t)));
attr = cell2struct({attr.val},strrep({attr.key},':','___'),2);
if ~all(ismember({'Version','NumberOfDataArrays'},fieldnames(attr)))
    error('[GIFTI] Missing mandatory attributes for GIFTI root element.');
end
if str2double(attr.Version) ~= 1
    warning('[GIFTI] Unknown specification version of GIFTI file (%s).',attr.Version);
end
nbData = str2double(attr.NumberOfDataArrays);

% Read children elements
%--------------------------------------------------------------------------
uid = children(t,root(t));
for i=1:length(uid)
    switch get(t,uid(i),'name')
        case 'MetaData'
            this.metadata    = gifti_MetaData(t,uid(i));
        case 'LabelTable'
            this.label       = gifti_LabelTable(t,uid(i));
        case 'DataArray'
            this.data{end+1} = gifti_DataArray(t,uid(i),filename);
        otherwise
            warning('[GIFTI] Unknown element "%s": ignored.',get(t,uid(i),'name'));
    end
end

if nbData ~= length(this.data)
    warning('[GIFTI] Mismatch between expected and effective number of datasets.');
end

%==========================================================================
function s = gifti_MetaData(t,uid)
s = struct('name',{}, 'value',{});
c = children(t,uid);
for i=1:length(c)
    for j=children(t,c(i))
        s(i).(lower(get(t,j,'name'))) = get(t,children(t,j),'value');
    end
end

%==========================================================================
function s = gifti_LabelTable(t,uid)
s = struct('name',{}, 'key',[], 'rgba',[]);
c = children(t,uid);
for i=1:length(c)
    a = attributes(t,'get',c(i));
    s(1).rgba(i,1:4) = NaN;
    for j=1:numel(a)
        switch lower(a{j}.key)
            case {'key','index'}
                s(1).key(i)    = str2double(a{j}.val);
            case 'red'
                s(1).rgba(i,1) = str2double(a{j}.val);
            case 'green'
                s(1).rgba(i,2) = str2double(a{j}.val);
            case 'blue'
                s(1).rgba(i,3) = str2double(a{j}.val);
            case 'alpha'
                s(1).rgba(i,4) = str2double(a{j}.val);
            otherwise
        end
    end
    s(1).name{i}  = get(t,children(t,c(i)),'value');
end

%==========================================================================
function s = gifti_DataArray(t,uid,filename)
s = struct(...
           'attributes', {}, ...
           'data',       {}, ...
           'metadata',   struct([]), ...
           'space',      {} ...
          );

attr = cell2mat(attributes(t,'get',uid));
s(1).attributes = cell2struct({attr.val},{attr.key},2);
s(1).attributes.Dim = [];
for i=1:str2double(s(1).attributes.Dimensionality)
    f = sprintf('Dim%d',i-1);
    s(1).attributes.Dim(i) = str2double(s(1).attributes.(f));
    s(1).attributes = rmfield(s(1).attributes,f);
end
s(1).attributes = rmfield(s(1).attributes,'Dimensionality');
if isfield(s(1).attributes,'ExternalFileName') && ...
        ~isempty(s(1).attributes.ExternalFileName)
    s(1).attributes.ExternalFileName = fullfile(fileparts(filename),...
        s(1).attributes.ExternalFileName);
end
    
c = children(t,uid);
for i=1:length(c)
    switch get(t,c(i),'name')
        case 'MetaData'
            s(1).metadata     = gifti_MetaData(t,c(i));
        case 'CoordinateSystemTransformMatrix'
            s(1).space(end+1) = gifti_Space(t,c(i));
        case 'Data'
            s(1).data         = gifti_Data(t,c(i),s(1).attributes);
        otherwise
            error('[GIFTI] Unknown DataArray element "%s".',get(t,c(i),'name'));
    end
end

%==========================================================================
function s = gifti_Space(t,uid)
s = struct('DataSpace','', 'TransformedSpace','', 'MatrixData',[]);
for i=children(t,uid)
    s.(get(t,i,'name')) = get(t,children(t,i),'value');
end
s.MatrixData = reshape(str2num(s.MatrixData),4,4)';

%==========================================================================
function d = gifti_Data(t,uid,s)
tp = getdict;
try
    tp = tp.(s.DataType);
catch
    error('[GIFTI] Unknown DataType.');
end
    
[unused,unused,mach] = fopen(1);
sb = @(x) x;
try
    if (strcmp(s.Endian,'LittleEndian') && ~isempty(strmatch('ieee-be',mach))) ...
        || (strcmp(s.Endian,'BigEndian') && ~isempty(strmatch('ieee-le',mach)))
        sb = @swapbyte;
    end
catch
    % Byte Order can be absent if encoding is ASCII, assume native otherwise
end

switch s.Encoding
    case 'ASCII'
        d = sscanf(get(t,children(t,uid),'value'),tp.format);

    case 'Base64Binary'
        d = typecast(sb(base64decode(get(t,children(t,uid),'value'))), tp.cast);

    case 'GZipBase64Binary'
        d = typecast(dunzip(sb(base64decode(get(t,children(t,uid),'value')))), tp.cast);

    case 'ExternalFileBinary'
        [p,f,e] = fileparts(s.ExternalFileName);
        if isempty(p)
            s.ExternalFileName = fullfile(pwd,[f e]);
        end
        if false
            fid = fopen(s.ExternalFileName,'r');
            if fid == -1
                error('[GIFTI] Unable to read binary file %s.',s.ExternalFileName);
            end
            fseek(fid,str2double(s.ExternalFileOffset),0);
            d = sb(fread(fid,prod(s.Dim),['*' tp.class]));
            fclose(fid);
        else
            d = file_array(s.ExternalFileName, s.Dim, tp.class, ...
                str2double(s.ExternalFileOffset),1,0,'rw');
        end
        
    otherwise
        error('[GIFTI] Unknown data encoding: %s.',s.Encoding);
end

if length(s.Dim) == 1, s.Dim(end+1) = 1; end
switch s.ArrayIndexingOrder
    case 'RowMajorOrder'
        d = permute(reshape(d,fliplr(s.Dim)),length(s.Dim):-1:1);
    case 'ColumnMajorOrder'
        d = reshape(d,s.Dim);
    otherwise
        error('[GIFTI] Unknown array indexing order.');
end
