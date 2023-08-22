function Header = spm_dicom_header(DicomFilename, DicomDictionary, Options)
% Read header information from a DICOM file
% FORMAT Header = spm_dicom_header(DicomFilename, DicomDictionary, Options)
% DicomFilename   - DICOM filename
% DicomDictionary - DICOM dictionary (see spm_dicom_headers)
% Options         - an (optional) structure containing fields
%                   abort      - if this is a function handle, it will
%                                be called with field name and value
%                                arguments.  If this function returns true,
%                                then reading the header will be aborted.
%                                [Default: false]
%                   all_fields - binary true/false, indicating what to do
%                                with fields that are not included in the
%                                DICOM dictionary.
%                                [Default: true]
%
% Header          - Contents of DICOM header
%
% Contents of headers are approximately explained in:
% http://medical.nema.org/standard.html
%
% This code may not work for all cases of DICOM data, as DICOM is an
% extremely complicated "standard".
%__________________________________________________________________________
% Copyright (C) 2002-2018 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dicom_header.m 7374 2018-07-09 17:09:46Z guillaume $


if nargin < 3
    Options = struct('abort',false, 'all_fields',true);
end

Header = [];
DicomFilename = deblank(DicomFilename);
FID  = fopen(DicomFilename,'r','ieee-le');
if FID == -1
    warning('spm:dicom','%s: Cant open file.', DicomFilename);
    return;
end

fseek(FID,128,'bof');
dcm = char(fread(FID,4,'uint8')');
if ~strcmp(dcm,'DICM')
    % Try truncated DICOM file fomat
    fseek(FID,0,'bof');
    Tag.Group   = fread(FID,1,'ushort');
    Tag.Element = fread(FID,1,'ushort');
    if isempty(Tag.Group) || isempty(Tag.Element)
        fclose(FID);
        warning('spm:dicom','%s: Truncated file.', DicomFilename);
        return;
    end
    if isempty(find(DicomDictionary.group==Tag.Group & DicomDictionary.element==Tag.Element,1)) && ~(Tag.Group==8 && Tag.Element==0)
        % Entry not found in DICOM dictionary and not from a GE Twin+excite
        % that starts with with an 8/0 Tag that I can't find any
        % documentation for.
        fclose(FID);
        warning('spm:dicom','%s: Not a DICOM file.', DicomFilename);
        return;
    else
        fseek(FID,0,'bof');
    end
end

try
    Header = ReadDicom(FID, 'il', DicomDictionary, Options);
    if isempty(Header)
        fclose(FID);
        return;
    end
    Header.Filename = fopen(FID);
catch
    le = lasterror;
    fprintf('%s: Trouble reading DICOM file (%s), skipping.\n', fopen(FID), le.message);
end

fclose(FID);


%==========================================================================
% function [Header, BytesRead] = ReadDicom(FID, TransferSyntax, DicomDictionary, Options, NumBytes)
%==========================================================================
function [Header, BytesRead] = ReadDicom(FID, TransferSyntax, DicomDictionary, Options, NumBytes)
if nargin<5, NumBytes = 4294967295; end % FFFFFFFF
BytesRead = 0;
Header    = [];
while BytesRead < NumBytes
    Tag = ReadTag(FID, TransferSyntax, DicomDictionary, Options);
    if isempty(Tag), break; end

    %fprintf('(%.4X,%.4X) "%s" %d %d %s\n', Tag.Group, Tag.Element, Tag.VR, Tag.Length, Tag.ExtraBytes, Tag.Name);

    if Tag.Group==65534 && Tag.Element==57357 % FFFE,E00D ItemDelimitationItem
        break;
    end
    if Tag.Length>0
        % Handle particular fields as special cases, otherwise just use default treatment
        switch Tag.Name
            case {'GroupLength'}
                % Ignore it
                fseek(FID,Tag.Length,'cof');
            case {'PixelData'}
                Header.StartOfPixelData = ftell(FID);
                Header.SizeOfPixelData  = Tag.Length;
                Header.VROfPixelData    = Tag.VR;
                fseek(FID,Tag.Length,'cof');
            case {'CSAData'} % raw data
                Header.StartOfCSAData = ftell(FID);
                Header.SizeOfCSAData  = Tag.Length;
                fseek(FID,Tag.Length,'cof');
            case {'CSAImageHeaderInfo', 'CSASeriesHeaderInfo','CSANonImageHeaderInfoVA',...
                  'CSAMiscProtocolHeaderInfoVA','CSANonImageHeaderInfoVB','CSAMiscProtocolHeaderInfoVB'}
                Content  = DecodeCSA(FID,Tag.Length);
                Header.(Tag.Name) = Content;
            case {'TransferSyntaxUID'}
                Content = deblank(char(fread(FID,Tag.Length,'uint8')'));
                Header.(Tag.Name) = Content;
                switch Content
                    case {'1.2.840.10008.1.2'}      % Implicit VR Little Endian
                        TransferSyntax = 'il';
                    case {'1.2.840.10008.1.2.1'}    % Explicit VR Little Endian
                        TransferSyntax = 'el';
                    case {'1.2.840.10008.1.2.1.99'} % Deflated Explicit VR Little Endian
                        warning('spm:dicom','%s: Cant read Deflated Explicit VR Little Endian file.', fopen(FID));
                       %TransferSyntax = 'dl';
                        return;
                    case {'1.2.840.10008.1.2.2'}    % Explicit VR Big Endian
                        %warning('spm:dicom','%s: Cant read Explicit VR Big Endian file',fopen(FID));
                        TransferSyntax = 'eb'; % Unused
                    case {'1.2.840.10008.1.2.4.50','1.2.840.10008.1.2.4.51','1.2.840.10008.1.2.4.70',...
                          '1.2.840.10008.1.2.4.80','1.2.840.10008.1.2.4.90','1.2.840.10008.1.2.4.91'} % JPEG Explicit VR
                        TransferSyntax = 'el';
                        %warning('spm:dicom',['Cant read JPEG Encoded file "' fopen(FID) '".']);
                    otherwise
                       %TransferSyntax = 'el';
                        warning('spm:dicom','%s: Unknown Transfer Syntax UID (%s).',fopen(FID), Content);
                        return;
                end
            otherwise
                % Default treatment
                switch Tag.VR
                    case {'UN'}
                        % Unknown - read as char
                        Content = fread(FID,Tag.Length,'uint8')';
                    case {'AE', 'AS', 'CS', 'DA', 'DS', 'DT', 'IS', 'LO', 'LT',...
                            'PN', 'SH', 'ST', 'TM', 'UI', 'UT'}
                        % Character strings
                        Content = char(fread(FID,Tag.Length,'uint8')');

                        switch Tag.VR
                            case {'UI','ST'}
                                Content = deblank(Content);
                            case {'DS'}
                                try
                                    Content = textscan(Content,'%f','delimiter','\\')';
                                    Content = Content{1};
                                catch
                                    Content = textscan(Content,'%f','delimiter','/')';
                                    Content = Content{1};
                                end
                            case {'IS'}
                                Content = textscan(Content,'%d','delimiter','\\')';
                                Content = double(Content{1});
                            case {'DA'}
                                Content = strrep(Content,'.',' ');
                                Content = textscan(Content,'%4d%2d%2d');
                                [y,m,d] = deal(Content{:});
                                Content = datenum(double(y),double(m),double(d));
                            case {'TM'}
                                if any(Content==':')
                                    Content = textscan(Content,'%d:%d:%f');
                                    [h,m,s] = deal(Content{:});
                                    h       = double(h);
                                    m       = double(m);
                                else
                                    Content = textscan(Content,'%2d%2d%f');
                                    [h,m,s] = deal(Content{:});
                                    h       = double(h);
                                    m       = double(m);
                                end
                                if isempty(h), h = 0; end
                                if isempty(m), m = 0; end
                                if isempty(s), s = 0; end
                                Content = s+60*(m+60*h);
                            case {'LO'}
                                Content = strrep(Content,'+AF8-','_');
                            otherwise
                        end
                    case {'OB'}
                        % dont know if this should be signed or unsigned
                        Content = fread(FID,Tag.Length,'uint8')';
                    case {'US', 'AT', 'OW'}
                        Content = fread(FID,Tag.Length/2,'uint16')';
                    case {'SS'}
                        Content = fread(FID,Tag.Length/2,'int16')';
                    case {'UL'}
                        Content = fread(FID,Tag.Length/4,'uint32')';
                    case {'SL'}
                        Content = fread(FID,Tag.Length/4,'int32')';
                    case {'FL','OF'}
                        Content = fread(FID,Tag.Length/4,'float')';
                    case {'FD','OD'}
                        Content = fread(FID,Tag.Length/8,'double')';
                    case {'SQ'}
                        [Content, BytesJustRead] = ReadSQ(FID, TransferSyntax, DicomDictionary, Options, Tag.Length);
                        if BytesJustRead==-1 && isempty(Content)
                            Header    = [];
                            BytesRead = -1;
                            return;
                        end
                        Tag.Length = BytesJustRead;
                    otherwise
                        Content = fread(FID,Tag.Length,'uint8')';
                end
                if ~isempty(Tag.Name)
                    if isa(Options.abort, 'function_handle') && feval(Options.abort, Tag.Name, Content)
                        Header    = [];
                        BytesRead = -1;
                        return;
                    end
                    Header.(Tag.Name) = Content;
                end
        end
    end
    BytesRead = BytesRead + Tag.ExtraBytes + Tag.Length;
end


%==========================================================================
% function [Header, BytesRead] = ReadSQ(FID, TransferSyntax, DicomDictionary, Options, NumBytes)
%==========================================================================
function [Header, BytesRead] = ReadSQ(FID, TransferSyntax, DicomDictionary, Options, NumBytes)
Header = {};
n      = 0;
BytesRead = 0;
while BytesRead<NumBytes
    Tag.Group   = fread(FID,1,'ushort');
    Tag.Element = fread(FID,1,'ushort');
    Tag.Length  = fread(FID,1,'uint');
    if isempty(Tag.Length), return; end % End of file
    BytesRead   = BytesRead + 8;
    if (Tag.Group == 65534) && (Tag.Element == 57344)   % FFFE/E000 Item
        [Item, BytesJustRead] = ReadDicom(FID, TransferSyntax, DicomDictionary, Options, Tag.Length);
        if BytesJustRead==-1 && isempty(Item)
            Header    = [];
            BytesRead = -1;
            return;
        end

        BytesRead    = BytesRead + BytesJustRead;
        if ~isempty(Item)
            n         = n + 1;
            Header{n} = Item;
        end
    elseif (Tag.Group == 65279) && (Tag.Element == 224) % FEFF/00E0 Item (Byte-swapped)
        % Byte-swapped
        [Filename,Permission,MachineFormat]      = fopen(FID);
        TransferSyntaxSwapped = TransferSyntax;
        if TransferSyntax(2)=='b'
            TransferSyntaxSwapped(2) = 'l';
        else
            TransferSyntaxSwapped(2) = 'b';
        end
        [Item, BytesJustRead] = ReadDicom(FID, TransferSyntaxSwapped, DicomDictionary, Options, Tag.Length);
        if BytesJustRead==-1 && isempty(Item)
            BytesRead = -1;
            return;
        end

        BytesRead = BytesRead + BytesJustRead;
        n         = n + 1;
        Header{n} = Item;
        Position  = ftell(FID);
        fclose(FID);
        FID     = fopen(Filename, Permission, MachineFormat);
        fseek(FID, Position, 'bof');
    elseif (Tag.Group == 65534) && (Tag.Element == 57565) % FFFE/E0DD SequenceDelimitationItem
        break;
    elseif (Tag.Group == 65279) && (Tag.Element == 56800) % FEFF/DDE0 SequenceDelimitationItem (Byte-swapped)
        % Byte-swapped
        break;
    else
        warning('spm:dicom','%s: Tag (%.4X,%.4X) is unexpected in sequence.', fopen(FID), Tag.Group, Tag.Element);
    end
end


%==========================================================================
% function Tag = ReadTag(FID,TransferSyntax,DicomDictionary, Options)
%==========================================================================
function Tag = ReadTag(FID, TransferSyntax, DicomDictionary, Options)
[GroupElement, Bytesread] = fread(FID,2,'ushort');
if Bytesread < 2, Tag = []; return; end
Tag.Group   = GroupElement(1);
Tag.Element = GroupElement(2);
if Tag.Group == 2, TransferSyntax = 'el'; end
t           = find(DicomDictionary.group==Tag.Group & DicomDictionary.element==Tag.Element);
if ~isempty(t)
    Tag.Name = DicomDictionary.values(t).name;
    Tag.VR   = DicomDictionary.values(t).vr{1};
else
    if ~Options.all_fields
        % With a reduced DicomDictionary, this can speed things up considerably.
        Tag.Name = '';
    else
        if Tag.Element~=0
            % Odd Group numbers indicate private (manufacturer-specific) fields
            if rem(Tag.Group,2)
                % Shadow Group
                Tag.Name = sprintf('Private_%.4x_%.4x',Tag.Group,Tag.Element);
            else
                Tag.Name = sprintf('Tag_%.4x_%.4x',Tag.Group,Tag.Element);
            end
            Tag.VR   = 'UN';
        else
            Tag.Name = '';
            Tag.VR   = 'UN';
        end
    end
end

if TransferSyntax(2) == 'b'
    [Filename,Permission,MachineFormat] = fopen(FID);
    if strcmp(MachineFormat,'ieee-le') || strcmp(MachineFormat,'ieee-le.l64')
        Position = ftell(FID);
        fclose(FID);
        FID      = fopen(Filename,Permission,'ieee-be');
        fseek(FID,Position,'bof');
    end
end

Tag.ExtraBytes = 4; % Two shorts read so far

if TransferSyntax(1) =='e'
    Tag.VR         = char(fread(FID,2,'uint8')');
    Tag.ExtraBytes = Tag.ExtraBytes + 2;
    switch Tag.VR
        case {'OB','OW','OF','OD','SQ','UN','UT'}
            if ~strcmp(Tag.VR,'UN') || Tag.Group~=65534
                fread(FID,1,'ushort');
                Tag.ExtraBytes = Tag.ExtraBytes + 2;
            else
                warning('spm:dicom','%s: Possible problem with %s Tag (VR="%s").', fopen(FID), Tag.Name, Tag.VR);
            end
            Tag.Length     = double(fread(FID,1,'uint'));
            Tag.ExtraBytes = Tag.ExtraBytes + 4;

        case {'AE','AS','AT','CS','DA','DS','DT','FD','FL','IS','LO','LT','PN','SH','SL','SS','ST','TM','UI','UL','US'}
            Tag.Length     = double(fread(FID,1,'ushort'));
            Tag.ExtraBytes = Tag.ExtraBytes + 2;

        case char([0 0])
            if (Tag.Group == 65534) && (Tag.Element == 57357)    % ItemDeliminationItem
                % at least on GE, ItemDeliminationItem does not have a
                % VR, but 4 bytes zeroes as Length
                Tag.VR         = 'UN';
                fread(FID,1,'ushort'); % Should be zero
                Tag.Length     = 0;
                Tag.ExtraBytes = Tag.ExtraBytes + 2;

            elseif (Tag.Group == 65534) && (Tag.Element == 57565) % SequenceDelimitationItem
                Tag.VR         = 'UN';
                Tag.Length     = 0;
                fread(FID,1,'ushort'); % Should be zero
                Tag.ExtraBytes = Tag.ExtraBytes + 2;

            else
                warning('spm:dicom','%s: Don''t know how to handle VR of "\0\0" in %s.',fopen(FID),Tag.Name);
            end
        otherwise
            warning('spm:dicom','%s: Possible problem with %s Tag (VR="%s")\n', fopen(FID), Tag.Name, Tag.VR);
            fread(FID,1,'ushort');
            Tag.Length     = double(fread(FID,1,'uint'));
            Tag.ExtraBytes = Tag.ExtraBytes + 2+4;
    end
else
    Tag.Length     = double(fread(FID,1,'uint'));
    Tag.ExtraBytes = Tag.ExtraBytes + 4;
end

if isempty(Tag.VR) || isempty(Tag.Length)
    Tag = [];
    return;
end


if rem(Tag.Length,2)
    if Tag.Length==4294967295 % FFFFFFFF
        if strcmp(Tag.VR,'UN')
            warning('spm:dicom','%s: Unknown Value Representation of %s Tag, assuming ''SQ''.', fopen(FID), Tag.Name);
            Tag.VR = 'SQ';
        end
        return;
    else
        warning('spm:dicom','%s: Odd numbered Value Length in %s Tag (%X).', fopen(FID), Tag.Name, Tag.Length);
        Tag = [];
    end
end



%==========================================================================
% function t = DecodeCSA(FID,NumBytes)
%==========================================================================
function t = DecodeCSA(FID,NumBytes)
% Decode shadow information (0029,1010) and (0029,1020) from Siemens files
[Filename,Permission,MachineFormat] = fopen(FID);
Position = ftell(FID);
if strcmp(MachineFormat,'ieee-be') || strcmp(MachineFormat,'ieee-be.l64')
    fclose(FID);
    FID  = fopen(Filename,Permission,'ieee-le');
    fseek(FID,Position,'bof');
end

c   = fread(FID,4,'uint8');
fseek(FID,Position,'bof');

if all(c'==[83 86 49 48]) % "SV10"
    t = DecodeCSA2(FID,NumBytes);
else
    t = DecodeCSA1(FID,NumBytes);
end

if strcmp(MachineFormat,'ieee-be') || strcmp(MachineFormat,'ieee-be.l64')
    fclose(FID);
    FID  = fopen(Filename, Permission, MachineFormat);
end
fseek(FID, Position+NumBytes, 'bof');


%==========================================================================
% function t = DecodeCSA1(FID,NumBytes)
%==========================================================================
function t = DecodeCSA1(FID,NumBytes)
n   = fread(FID,1,'uint32');
if isempty(n) || n>1024 || n <= 0
    fseek(FID,NumBytes-4,'cof');
    t = struct('name','JUNK: Don''t know how to read this damned file format');
    return;
end
fread(FID,1,'uint32'); % Unused "M" or 77 for some reason
BytesRead = 2*4;
t(n)      = struct('name','', 'vm',[], 'vr','', 'syngodt',[], 'nitems',[], 'xx',[], 'item',struct('xx',{},'val',{}));
for i=1:n
    t(i).name = fread(FID,64,'uint8')';
    msk       = find(~t(i).name)-1;
    if ~isempty(msk)
        t(i).name = char(t(i).name(1:msk(1)));
    else
        t(i).name = char(t(i).name);
    end
    t(i).vm      = fread(FID,1,'int32')';
    t(i).vr      = fread(FID,4,'uint8')';
    t(i).vr      = char(t(i).vr(1:3));
    t(i).syngodt = fread(FID,1,'int32')';
    t(i).nitems  = fread(FID,1,'int32')';
    t(i).xx      = fread(FID,1,'int32')'; % 77 or 205
    BytesRead    = BytesRead + 64+4+4+4+4+4;
    if t(i).nitems > 0
        t(i).item(t(i).nitems) = struct('xx',[], 'val',[]);
    end
    for j=1:t(i).nitems
        % This bit is just wierd
        t(i).item(j).xx  = fread(FID,4,'int32')'; % [x x 77 x]
        BytesToRead      = t(i).item(j).xx(1)-t(1).nitems;
        if BytesToRead<0 || BytesToRead+Bytesread+4*4>NumBytes
            t(i).item(j).val = '';
            t(i).item        = t(i).item(1:j);
            BytesRead        = BytesRead + 4*4;
            break;
        end
        t(i).item(j).val = char(fread(FID,BytesToRead,'uint8')');
        fread(FID,4-rem(BytesToRead,4),'uint8');
        BytesRead        = BytesRead + 4*4+BytesToRead+(4-rem(BytesToRead,4));
    end
end


%==========================================================================
% function t = DecodeCSA2(FID, NumBytes)
%==========================================================================
function t = DecodeCSA2(FID, NumBytes)
fread(FID, 8, 'uint8'); % Unused
n        = fread(FID, 1, 'uint32');
if isempty(n) || n>1024 || n < 0
    fseek(FID, NumBytes-4, 'cof');
    t = struct('name','Don''t know how to read this damned file format');
    return;
end
fread(FID, 1, 'uint32'); % Unused "M" or 77 for some reason
Position = 16;
t(n)     = struct('name','', 'vm',[], 'vr','', 'syngodt',[], 'nitems',[], 'xx',[], 'item',struct('xx',{},'val',{}));
for i=1:n
    t(i).name = fread(FID, 64, 'uint8')';
    Position  = Position + 64;
    msk       = find(~t(i).name)-1;
    if ~isempty(msk)
        t(i).name = char(t(i).name(1:msk(1)));
    else
        t(i).name = char(t(i).name);
    end
    t(i).vm      = fread(FID, 1, 'int32')';
    t(i).vr      = fread(FID, 4, 'uint8')';
    t(i).vr      = char(t(i).vr(1:3));
    t(i).syngodt = fread(FID, 1, 'int32')';
    t(i).nitems  = fread(FID, 1, 'int32')';
    t(i).xx      = fread(FID, 1, 'int32')'; % 77 or 205
    Position     = Position + 20;
    if t(i).nitems > 0
        t(i).item(t(i).nitems) = struct('xx',[], 'val',[]);
    end
    for j=1:t(i).nitems
        t(i).item(j).xx  = fread(FID,4,'int32')'; % [x x 77 x]
        Position         = Position + 16;
        BytesToRead      = t(i).item(j).xx(2);
        if BytesToRead > NumBytes-Position
            BytesToRead      = NumBytes-Position;
            t(i).item(j).val = char(fread(FID, BytesToRead, 'uint8')');
            fread(FID,rem(4-rem(BytesToRead, 4), 4), 'uint8');
            t(i).item        = t(i).item(1:j);
            warning('spm:dicom','%s: Problem reading Siemens CSA field.', fopen(FID));
            return;
        end
        t(i).item(j).val = char(fread(FID, BytesToRead, 'uint8')');
        fread(FID,rem(4-rem(BytesToRead, 4), 4), 'uint8');
    end
end
