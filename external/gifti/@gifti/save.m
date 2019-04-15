function save(this,filename,encoding)
% Save GIfTI object in a GIfTI format file
% FORMAT save(this,filename,encoding)
% this      - GIfTI object
% filename  - name of GIfTI file to be created [Default: 'untitled.gii']
% encoding  - optional argument to specify encoding format, among
%             ASCII, Base64Binary, GZipBase64Binary, ExternalFileBinary.
%             [Default: 'GZipBase64Binary']
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: save.m 6516 2015-08-07 17:28:33Z guillaume $


% Check filename
%--------------------------------------------------------------------------
ext = '.gii';
if nargin == 1
    filename = 'untitled';
end
[p,f,e] = fileparts(filename);
if ~ismember(lower(e),{ext})
    e = ext;
end
filename = fullfile(p,[f e]);

% Open file for writing
%--------------------------------------------------------------------------
fid = fopen(filename,'wt');
if fid == -1
    error('Unable to write file %s: permission denied.',filename);
end

% Write file
%--------------------------------------------------------------------------
if nargin < 3, encoding = 'GZipBase64Binary'; end
switch encoding
    case {'ASCII','Base64Binary','GZipBase64Binary','ExternalFileBinary'}
    otherwise
        error('Unknown encoding.');
end
fid = save_gii(fid,this,encoding);

% Close file
%--------------------------------------------------------------------------
fclose(fid);


%==========================================================================
% function fid = save_gii(fid,this,encoding)
%==========================================================================
function fid = save_gii(fid,this,encoding)

% Defaults for DataArray's attributes
%--------------------------------------------------------------------------
[unused,unused,mach]   = fopen(fid);
if strncmp('ieee-be',mach,7)
    def.Endian         = 'BigEndian';
elseif strncmp('ieee-le',mach,7)
    def.Endian         = 'LittleEndian';
else
    error('[GIFTI] Unknown byte order "%s".',mach);
end
def.Encoding           = encoding;
def.Intent             = 'NIFTI_INTENT_NONE';
def.DataType           = 'NIFTI_TYPE_FLOAT32';
def.ExternalFileName   = '';
def.ExternalFileOffset = '';
def.offset             = 0;

% Edit object DataArray attributes
%--------------------------------------------------------------------------
for i=1:length(this.data)
    % Revert the dimension storage
    d = this.data{i}.attributes.Dim;
    if numel(d) > 1 && d(end) == 1
        d = d(1:end-1);
    end
    this.data{i}.attributes = rmfield(this.data{i}.attributes,'Dim');
    this.data{i}.attributes.Dimensionality = num2str(length(d));
    for j=1:length(d)
        this.data{i}.attributes.(sprintf('Dim%d',j-1)) = num2str(d(j));
    end
    % Enforce some conventions
    this.data{i}.attributes.ArrayIndexingOrder = 'ColumnMajorOrder';
    if ~isfield(this.data{i}.attributes,'DataType') || ...
        isempty(this.data{i}.attributes.DataType)
        warning('DataType set to default: %s', def.DataType);
        this.data{i}.attributes.DataType       = def.DataType;
    end
    if ~isfield(this.data{i}.attributes,'Intent') || ...
        isempty(this.data{i}.attributes.Intent)
        warning('Intent code set to default: %s', def.Intent);
        this.data{i}.attributes.Intent         = def.Intent;
    end
    this.data{i}.attributes.Encoding           = def.Encoding;
    this.data{i}.attributes.Endian             = def.Endian;
    this.data{i}.attributes.ExternalFileName   = def.ExternalFileName;
    this.data{i}.attributes.ExternalFileOffset = def.ExternalFileOffset;
    switch this.data{i}.attributes.Encoding
        case {'ASCII', 'Base64Binary','GZipBase64Binary' }
        case 'ExternalFileBinary'
            extfilename = this.data{i}.attributes.ExternalFileName;
            if isempty(extfilename)
                [p,f] = fileparts(fopen(fid));
                extfilename = [f '.dat'];
            end
            [p,f,e] = fileparts(extfilename);
            this.data{i}.attributes.ExternalFileName   = fullfile(fileparts(fopen(fid)),[f e]);
            this.data{i}.attributes.ExternalFileOffset = num2str(def.offset);
        otherwise
            error('[GIFTI] Unknown data encoding: %s.',this.data{i}.attributes.Encoding);
    end
end

% Prolog
%--------------------------------------------------------------------------
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<!DOCTYPE GIFTI SYSTEM "http://www.nitrc.org/frs/download.php/115/gifti.dtd">\n');
fprintf(fid,'<GIFTI Version="1.0"  NumberOfDataArrays="%d">\n',numel(this.data));

o = @(x) blanks(x*3);

% MetaData
%--------------------------------------------------------------------------
fprintf(fid,'%s<MetaData',o(1));
if isempty(this.metadata)
    fprintf(fid,'/>\n');
else
    fprintf(fid,'>\n');
    for i=1:length(this.metadata)
        fprintf(fid,'%s<MD>\n',o(2));
        fprintf(fid,'%s<Name><![CDATA[%s]]></Name>\n',o(3),...
            this.metadata(i).name);
        fprintf(fid,'%s<Value><![CDATA[%s]]></Value>\n',o(3),...
            this.metadata(i).value);
        fprintf(fid,'%s</MD>\n',o(2));
    end
    fprintf(fid,'%s</MetaData>\n',o(1));
end

% LabelTable
%--------------------------------------------------------------------------
fprintf(fid,'%s<LabelTable',o(1));
if isempty(this.label)
    fprintf(fid,'/>\n');
else
    fprintf(fid,'>\n');
    for i=1:length(this.label.name)
        if ~all(isnan(this.label.rgba(i,:)))
            label_rgba = sprintf(' Red="%f" Green="%f" Blue="%f" Alpha="%f"',...
                this.label.rgba(i,:));
        else
            label_rgba = '';
        end
        fprintf(fid,'%s<Label Key="%d"%s><![CDATA[%s]]></Label>\n',o(2),...
            this.label.key(i), label_rgba, this.label.name{i});
    end
    fprintf(fid,'%s</LabelTable>\n',o(1));
end

% DataArray
%--------------------------------------------------------------------------
for i=1:length(this.data)
    fprintf(fid,'%s<DataArray',o(1));
    if def.offset
        this.data{i}.attributes.ExternalFileOffset = num2str(def.offset);
    end
    fn = sort(fieldnames(this.data{i}.attributes));
    oo = repmat({o(5) '\n'},length(fn),1); oo{1} = '  '; oo{end} = '';
    for j=1:length(fn)
        if strcmp(fn{j},'ExternalFileName')
            [p,f,e] = fileparts(this.data{i}.attributes.(fn{j}));
            attval = [f e];
        else
            attval = this.data{i}.attributes.(fn{j});
        end
        fprintf(fid,'%s%s="%s"%s',oo{j,1},...
                fn{j},attval,sprintf(oo{j,2}));
    end
    fprintf(fid,'>\n');
    
    % MetaData
    %----------------------------------------------------------------------
    fprintf(fid,'%s<MetaData>\n',o(2));
    for j=1:length(this.data{i}.metadata)
        fprintf(fid,'%s<MD>\n',o(3));
        fprintf(fid,'%s<Name><![CDATA[%s]]></Name>\n',o(4),...
            this.data{i}.metadata(j).name);
        fprintf(fid,'%s<Value><![CDATA[%s]]></Value>\n',o(4),...
            this.data{i}.metadata(j).value);
        fprintf(fid,'%s</MD>\n',o(3));
    end
    fprintf(fid,'%s</MetaData>\n',o(2));
    
    % CoordinateSystemTransformMatrix
    %----------------------------------------------------------------------
    for j=1:length(this.data{i}.space)
        fprintf(fid,'%s<CoordinateSystemTransformMatrix>\n',o(2));
        fprintf(fid,'%s<DataSpace><![CDATA[%s]]></DataSpace>\n',o(3),...
            this.data{i}.space(j).DataSpace);
        fprintf(fid,'%s<TransformedSpace><![CDATA[%s]]></TransformedSpace>\n',o(3),...
            this.data{i}.space(j).TransformedSpace);
        fprintf(fid,'%s<MatrixData>%s</MatrixData>\n',o(3),...
            sprintf('%f ',this.data{i}.space(j).MatrixData'));
        fprintf(fid,'%s</CoordinateSystemTransformMatrix>\n',o(2));
    end
    
    % Data (saved using MATLAB's ColumnMajorOrder)
    %----------------------------------------------------------------------
    fprintf(fid,'%s<Data>',o(2));
    tp = getdict;
    try
        tp = tp.(this.data{i}.attributes.DataType);
    catch
        error('[GIFTI] Unknown DataType.');
    end
    switch this.data{i}.attributes.Encoding
        case 'ASCII'
            fprintf(fid, [tp.format ' '], this.data{i}.data);
        case 'Base64Binary'
            fprintf(fid,base64encode(typecast(this.data{i}.data(:),'uint8')));
            % uses native machine format
        case 'GZipBase64Binary'
            fprintf(fid,base64encode(zstream('C',typecast(this.data{i}.data(:),'uint8'))));
            % uses native machine format
        case 'ExternalFileBinary'
            extfilename = this.data{i}.attributes.ExternalFileName;
            dat = this.data{i}.data;
            if isa(dat,'file_array')
                dat = subsref(dat,substruct('()',repmat({':'},1,numel(dat.dim))));
            end
            if ~def.offset
                fide = fopen(extfilename,'w'); % uses native machine format
            else
                fide = fopen(extfilename,'a'); % uses native machine format
            end
            if fide == -1
                error('Unable to write file %s: permission denied.',extfilename);
            end
            fseek(fide,0,1);
            fwrite(fide,dat,tp.class);
            def.offset = ftell(fide);
            fclose(fide);
        otherwise
            error('[GIFTI] Unknown data encoding.');
    end
    fprintf(fid,'</Data>\n');
    fprintf(fid,'%s</DataArray>\n',o(1));
end

fprintf(fid,'</GIFTI>\n');
