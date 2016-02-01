function [data] = fiff_setup_read_raw(fname,allow_maxshield)
%
% [data] = fiff_setup_read_raw(fname,allow_maxshield)
%
% Read information about raw data file
%
% fname               Name of the file to read
% allow_maxshield     Accept unprocessed MaxShield data
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.12  2009/03/30 11:37:37  msh
%   Added copying of measurement info blocks from the original like in mne_browse_raw
%
%   Revision 1.11  2008/09/30 17:22:39  msh
%   Added possibility to load unprocessed MaxShield data
%
%   Revision 1.10  2008/06/16 20:26:03  msh
%   Added some semicolons to avoid unncessary output.
%
%   Revision 1.9  2008/06/16 19:48:32  msh
%   Added some comments
%
%   Revision 1.8  2008/06/16 19:46:30  msh
%   Fixed handling of initial skip so that the result is compatible with the event files
%
%   Revision 1.7  2007/06/08 14:44:16  msh
%   Fixed problem appearing with multiple skips in a data file
%
%   Revision 1.6  2006/05/03 18:53:05  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.5  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.4  2006/04/21 16:17:48  msh
%   Added handling of CTF compensation
%
%   Revision 1.3  2006/04/21 14:23:16  msh
%   Further improvements in raw data reading
%
%   Revision 1.2  2006/04/21 11:33:18  msh
%   Added raw segment reading routines
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_setup_read_raw';

if nargin ~= 1 && nargin ~= 2
    error(me,'Incorrect number of arguments');
end

if nargin == 1
    allow_maxshield = false;
end
%
%   Open the file
%
fprintf(1,'Opening raw data file %s...\n',fname);
[ fid, tree ] = fiff_open(fname);
%
%   Read the measurement info
%
[ info, meas ] = fiff_read_meas_info(fid,tree);
%
%   Locate the data of interest
%
raw = fiff_dir_tree_find(meas,FIFF.FIFFB_RAW_DATA);
if isempty(raw)
    raw = fiff_dir_tree_find(meas,FIFF.FIFFB_CONTINUOUS_DATA);
    if allow_maxshield
        raw = fiff_dir_tree_find(meas,FIFF.FIFFB_SMSH_RAW_DATA);
        if isempty(raw)
            error(me,'No raw data in %s',fname);
        end
    else
        if isempty(raw)
            error(me,'No raw data in %s',fname);
        end
    end
end
%
%   Set up the output structure
%
info.filename   = fname;
data.fid        = fid;
data.info       = info;
data.first_samp = 0;
data.last_samp  = 0;
%
%   Process the directory
%
dir          = raw.dir;
nent         = raw.nent;
nchan        = info.nchan;
first        = 1;
first_samp   = 0;
first_skip   = 0;
%
%  Get first sample tag if it is there
%
if dir(first).kind == FIFF.FIFF_FIRST_SAMPLE
    tag = fiff_read_tag(fid,dir(first).pos);
    first_samp = tag.data;
    first = first + 1;
end
%
%  Omit initial skip
%
if dir(first).kind == FIFF.FIFF_DATA_SKIP
    %
    %  This first skip can be applied only after we know the buffer size
    %
    tag = fiff_read_tag(fid,dir(first).pos);
    first_skip = tag.data;
    first = first + 1;
end
%
%  Get first sample tag if it is there
%
if dir(first).kind == FIFF.FIFF_FIRST_SAMPLE
    tag = fiff_read_tag(fid,dir(first).pos);
    first_samp = first_samp + tag.data;
    first = first + 1;
end
data.first_samp = first_samp;
%
%   Go through the remaining tags in the directory
%
rawdir = struct('ent',{},'first',{},'last',{},'nsamp',{});
nskip = 0;
ndir  = 0;
for k = first:nent
    ent = dir(k);
    if ent.kind == FIFF.FIFF_DATA_SKIP
        tag = fiff_read_tag(fid,ent.pos);
        nskip = tag.data;
    elseif ent.kind == FIFF.FIFF_DATA_BUFFER
        %
        %   Figure out the number of samples in this buffer
        %
        switch ent.type
            case FIFF.FIFFT_DAU_PACK16
                nsamp = ent.size/(2*nchan);
            case FIFF.FIFFT_SHORT
                nsamp = ent.size/(2*nchan);
            case FIFF.FIFFT_FLOAT
                nsamp = ent.size/(4*nchan);
            case FIFF.FIFFT_INT
                nsamp = ent.size/(4*nchan);
            otherwise
                fclose(fid);
                error(me,'Cannot handle data buffers of type %d',ent.type);
        end
        %
        %  Do we have an initial skip pending?
        %
        if first_skip > 0
            first_samp = first_samp + nsamp*first_skip;
            data.first_samp = first_samp;
            first_skip = 0;
        end
        %
        %  Do we have a skip pending?
        %
        if nskip > 0
            ndir        = ndir+1;
            rawdir(ndir).ent   = [];
            rawdir(ndir).first = first_samp;
            rawdir(ndir).last  = first_samp + nskip*nsamp - 1;
            rawdir(ndir).nsamp = nskip*nsamp;
            first_samp = first_samp + nskip*nsamp;
            nskip = 0;
        end
        %
        %  Add a data buffer
        %
        ndir               = ndir+1;
        rawdir(ndir).ent   = ent;
        rawdir(ndir).first = first_samp;
        rawdir(ndir).last  = first_samp + nsamp - 1;
        rawdir(ndir).nsamp = nsamp;
        first_samp = first_samp + nsamp;
    end
end
data.last_samp  = first_samp - 1;
%
%   Add the calibration factors
%
cals = zeros(1,data.info.nchan);
for k = 1:data.info.nchan
    cals(k) = data.info.chs(k).range*data.info.chs(k).cal;
end
%
data.cals       = cals;
data.rawdir     = rawdir;
data.proj       = [];
data.comp       = [];
%
fprintf(1,'\tRange : %d ... %d  =  %9.3f ... %9.3f secs\n',...
    data.first_samp,data.last_samp,...
    double(data.first_samp)/data.info.sfreq,...
    double(data.last_samp)/data.info.sfreq);
fprintf(1,'Ready.\n');
fclose(data.fid);
data.fid = -1;
return;

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
