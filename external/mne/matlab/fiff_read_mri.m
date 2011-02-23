function [stack] = fiff_read_mri(fname,read_data)
%
% [stack] = fiff_read_mri(fname,read_data)
%
% read_data argument is optional, if set to false the pixel data are
% not read. The default is to read the pixel data
%
% Read a fif format MRI description file
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.5  2009/04/05 13:56:28  msh
%   Fixed a typo in the log messages
%
%   Revision 1.4  2009/04/05 13:50:08  msh
%   Added voxel -> MRI coordinate transform
%
%   Revision 1.3  2009/01/18 23:50:52  msh
%   Handle ushort embedded data correctly
%
%   Revision 1.2  2008/11/17 00:23:51  msh
%   Fixed error in the coordinate system transformation between voxel and RAS coordinates
%
%   Revision 1.1  2006/04/26 19:50:58  msh
%   Added fiff_read_mri
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_read_mri';

if nargin == 1
    read_data = true;
elseif nargin ~= 2
    error(me,'Incorrect number of arguments');
end
%
%   Try to open the file
%
[ fid, tree ] = fiff_open(fname);
%
%   Locate the data of interest
%   Pick the first MRI set within the first MRI data block
%
mri = fiff_dir_tree_find(tree,FIFF.FIFFB_MRI);
if length(mri) == 0
    fclose(fid);
    error(me,'Could not find MRI data');
end
mri = mri(1);
%
set = fiff_dir_tree_find(mri,FIFF.FIFFB_MRI_SET);
if length(set) == 0
    fclose(fid);
    error(me,'Could not find MRI stack');
end
set = set(1);
%
slices = fiff_dir_tree_find(set,FIFF.FIFFB_MRI_SLICE);
if length(slices) == 0
    fclose(fid);
    error(me,'Could not find MRI slices');
end
%
%   Data ID
%
tag = find_tag(mri,FIFF.FIFF_BLOCK_ID);
if ~isempty(tag)
    stack.id = tag.data;
end
%
%   Head -> MRI coordinate transformation
%   MRI (surface RAS) -> MRI (RAS) transformation
%
stack.trans = [];
stack.ras_trans = [];
for p = 1:set.nent
    kind = set.dir(p).kind;
    if kind == FIFF.FIFF_COORD_TRANS
        tag = fiff_read_tag(fid,set.dir(p).pos);
        if tag.data.to == FIFF.FIFFV_COORD_MRI && tag.data.from == FIFF.FIFFV_COORD_HEAD
            stack.trans = tag.data;
        elseif tag.data.to == FIFF.FIFFV_MNE_COORD_RAS && tag.data.from == FIFF.FIFFV_COORD_MRI
            stack.ras_trans = tag.data;
        end
    end
end
if isempty(stack.trans)
    stack.trans.from  = FIFF.FIFFV_COORD_HEAD;
    stack.trans.to    = FIFF.FIFFV_COORD_MRI;
    stack.trans.trans = eye(4,4);
end
if isempty(stack.ras_trans)
    stack.ras_trans.from   = FIFF.FIFFV_COORD_MRI;
    stack.ras_trans.to     = FIFF.FIFFV_MNE_COORD_RAS;
    stack.ras_trans.trans = eye(4,4);
end
stack.voxel_trans = [];
stack.nslice = length(slices);
if read_data
    fprintf(1,'\tReading slice information and pixel data.');
else
    fprintf(1,'\tReading slice information.');
end
for k = 1:stack.nslice
    try
        stack.slices(k) = read_slice(slices(k));
    catch
        fclose(fid);
        error(me,'%s',mne_omit_first_line(lasterr));
    end
    if mod(k,50) == 0
        fprintf(1,'.%d.',k);
    end
end
fprintf(1,'.%d..[done]\n',k);

add_voxel_transform

fclose(fid);


return;

    function add_voxel_transform
        
        if stack.nslice < 2
            fprintf(1,'\tOnly one slice in the stack. Voxel transformation will not be included\n');
            return;
        end
        %
        %   Check that the slices really form a stack
        %
        for k = 1:stack.nslice-1
            d = stack.slices(k+1).trans.trans(1:3,4)-stack.slices(k).trans.trans(1:3,4);
            d = sqrt(d'*d);
            if k == 1
                d0 = d;
                r0 = stack.slices(k).trans.trans(1:3,4);
                n0 = stack.slices(k).trans.trans(1:3,3);
            else
                if abs(stack.slices(k).pixel_width-stack.slices(1).pixel_width) > 1e-4 || abs(stack.slices(k).pixel_height-stack.slices(1).pixel_height) > 1e-4
                    fprintf(1,'\tPixel sizes are not equal. Voxel transformation will not be included\n');
                    return;
                end
                if stack.slices(k).width ~= stack.slices(1).width || stack.slices(k).height ~= stack.slices(1).height
                    fprintf(1,'\tImage sizes are not equal. Voxel transformation will not be included\n');
                    return;
                end
                if abs(d-d0) > 1e-4
                    fprintf(1,'\tThe slices are not equally spaced. Voxel transformation will not be included\n');
                    return;
                end
                %
                %   Rectangular volume?
                %
                r1 = r0 + (k-1)*d0*n0;
                d = stack.slices(k).trans.trans(1:3,4) - r1;
                d = sqrt(d'*d);
                if abs(d) > 1e-4
                    fprintf(1,'\tThe slices do not form a rectangular volume. Voxel transformation will not be included\n');
                    return;
                end
            end
        end
        %
        %   Ready to proceed
        %
        voxel_trans = stack.slices(1).trans;
        t = eye(4);
        t(1,1) = stack.slices(1).pixel_width;
        t(2,2) = stack.slices(1).pixel_height;
        d = stack.slices(2).trans.trans(1:3,4)-stack.slices(1).trans.trans(1:3,4);
        t(3,3) = sqrt(d'*d);
        voxel_trans.trans = voxel_trans.trans*t;
        voxel_trans.from  = FIFF.FIFFV_MNE_COORD_MRI_VOXEL;
        stack.voxel_trans = voxel_trans;
        
        fprintf(1,'\tVoxel transformation added\n');
        
    end

    function [slice] = read_slice(node)
        %
        %   Read all components of a single slice
        %
        tag = find_tag(node,FIFF.FIFF_COORD_TRANS);
        if isempty(tag)
            error(me,'Could not find slice coordinate transformation');
        end
        slice.trans = tag.data;
        if slice.trans.from ~= FIFF.FIFFV_COORD_MRI_SLICE || ...
                slice.trans.to ~= FIFF.FIFFV_COORD_MRI
            error(me,'Illegal slice coordinate transformation');
        end
        %
        %   Change the coordinate transformation so that
        %   ex is right
        %   ey is down (up in the fif file)
        %   ez steps to the next slice in the series
        %
        slice.trans.trans(1:3,2) = -slice.trans.trans(1:3,2);
        %
        %   Pixel data info
        %
        tag = find_tag(node,FIFF.FIFF_MRI_PIXEL_ENCODING);
        if isempty(tag)
            error(me,'Pixel encoding tag missing');
        end
        slice.encoding = tag.data;
        %
        %    Offset in the MRI data file if not embedded
        %
        tag = find_tag(node,FIFF.FIFF_MRI_PIXEL_DATA_OFFSET);
        if isempty(tag)
            slice.offset = -1;
        else
            slice.offset = tag.data;
        end
        %
        %    Format of the MRI source file
        %
        tag = find_tag(node,FIFF.FIFF_MRI_SOURCE_FORMAT);
        if isempty(tag)
            slice.source_format = 0;
        else
            slice.source_format = tag.data;
        end
        %
        %   Suggested scaling for the pixel values
        %   (not applied here)
        %
        tag = find_tag(node,FIFF.FIFF_MRI_PIXEL_SCALE);
        if isempty(tag)
            slice.scale = 1.0;
        else
            slice.scale = tag.data;
        end
        %
        %   Width and height in pixels
        %
        tag = find_tag(node,FIFF.FIFF_MRI_WIDTH);
        if isempty(tag)
            error(me,'Slice width missing');
        end
        slice.width = tag.data;
        %
        tag = find_tag(node,FIFF.FIFF_MRI_HEIGHT);
        if isempty(tag)
            error(me,'Slice height missing');
        end
        slice.height = tag.data;
        %
        %   Pixel sizes
        %
        tag = find_tag(node,FIFF.FIFF_MRI_WIDTH_M);
        if isempty(tag)
            error(me,'Pixel width missing');
        end
        slice.pixel_width = double(tag.data)/double(slice.width);
        %
        tag = find_tag(node,FIFF.FIFF_MRI_HEIGHT_M);
        if isempty(tag)
            error(me,'Pixel height missing');
        end
        slice.pixel_height = double(tag.data)/double(slice.height);
        %
        %   Are the data here or in another file?
        %
        tag = find_tag(node,FIFF.FIFF_MRI_SOURCE_PATH);
        if isempty(tag)
            slice.offset = -1;
            slice.source = [];
            %
            %   Pixel data are embedded in the fif file
            %
            if read_data
                tag = find_tag(node,FIFF.FIFF_MRI_PIXEL_DATA);
                if isempty(tag)
                    error(me,'Embedded pixel data missing');
                end
                if slice.encoding == FIFF.FIFFV_MRI_PIXEL_WORD
                    if tag.type ~= slice.encoding && tag.type ~= FIFF.FIFFT_USHORT
                        error(me,'Embedded data is in wrong format (expected %d, got %d)',...
                            slice.encoding,tag.type);
                    end
                else
                    if tag.type ~= slice.encoding
                        error(me,'Embedded data is in wrong format (expected %d, got %d)',...
                            slice.encoding,tag.type);
                    end
                end
                if length(tag.data) ~= slice.width*slice.height
                    error(me,'Wrong length of pixel data');
                end
                %
                %   Reshape into an image
                %
                slice.data = reshape(tag.data,slice.width,slice.height)';
            end
        else
            if slice.offset < 0
                error(me,'Offset to external file missing');
            end
            slice.source = tag.data;
            %
            %   External slice reading follows
            %
            if read_data
                pname = search_pixel_file(slice.source,fname);
                if isempty(pname)
                    error(me,'Could not locate pixel file %s',slice.source);
                else
                    try
                        slice.data   = read_external_pixels(pname,...
                            slice.offset,slice.encoding,slice.width,slice.height);
                    catch
                        error(me,'%s',mne_omit_first_line(lasterr));
                    end
                end
            end
        end
    end

    function [name] = search_pixel_file(pname,sname)
        %
        %   First try the file name as it is
        %
        if exist(pname,'file') == 2
            name = pname;
        else
            %
            %   Then <set file dir>/../slices/<slice file name>
            %
            a = findstr(sname,'/');
            if isempty(a)
                d = pwd;
            else
                d = sname(1:a(length(a))-1);
            end
            a = findstr(pname,'/');
            if ~isempty(a)
                pname = pname(a(length(a))+1:length(pname));
            end
            pname = sprintf('%s/../slices/%s',d,pname);
            if exist(pname,'file') == 2
                name = pname;
            else
                name = [];
            end
        end
        return;
        
    end

    function [pixels] = read_external_pixels(pname,offset,encoding,width,height)
        %
        %   Read pixel data from an external file
        %
        if (encoding == FIFF.FIFFV_MRI_PIXEL_SWAP_WORD)
            sfid = fopen(pname,'rb','ieee-le');
        else
            sfid = fopen(pname,'rb','ieee-be');
        end
        if sfid < 0
            error(me,'Could not open pixel data file : %s',pname);
        end
        try
            fseek(sfid,double(offset),'bof');
        catch
            error(me,'Could not position to pixel data @ %d',offset);
            fclose(sfid);
        end
        %
        %   Proceed carefully according to the encoding
        %
        switch encoding
            
            case FIFF.FIFFV_MRI_PIXEL_BYTE
                pixels = fread(sfid,double(width*height),'uint8=>uint8');
            case FIFF.FIFFV_MRI_PIXEL_WORD
                pixels = fread(sfid,double(width*height),'int16=>int16');
            case FIFF.FIFFV_MRI_PIXEL_SWAP_WORD
                pixels = fread(sfid,double(width*height),'int16=>int16');
            case FIFF.FIFFV_MRI_PIXEL_FLOAT
                pixels = fread(sfid,double(width*height),'single=>double');
            otherwise
                fclose(sfid);
                error(me,'Unknown pixel encoding : %d',encoding);
        end
        fclose(sfid);
        %
        %   Reshape into an image
        %
        pixels = reshape(pixels,width,height)';
    end

    function [tag] = find_tag(node,findkind)
        
        for p = 1:node.nent
            kind = node.dir(p).kind;
            pos  = node.dir(p).pos;
            if kind == findkind
                tag = fiff_read_tag(fid,pos);
                return;
            end
        end
        tag = [];
        return
        
    end

end

