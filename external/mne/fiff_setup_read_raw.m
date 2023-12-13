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
end
%
%   Special handling of data recorded with Internal Active Shielding
%
if isempty(raw) && allow_maxshield
    raw = fiff_dir_tree_find(meas,FIFF.FIFFB_IAS_RAW_DATA);
    if ~isempty(raw)
        disp([10 '--------' 10 ...
        'WARNING: This file contains raw Internal Active Shielding data. It may be distorted.' 10 ...
        'Elekta recommends it be run through MaxFilter to produce reliable results.' 10 ...
        'Consider closing the file and running MaxFilter on the data.' 10 ...
        '--------' 10]);
    end
end
if isempty(raw)
    error(me,'No raw data in file');
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
%   Check whether the data are organized uniformly as one sample per tag,
%   which may create an administrative overhead in file handling, that can
%   be avoided. The extraction of the data from the struct array may take
%   some time for large arrays, but is worth the wait, because still faster
%   than looping across all tags in the other case. Note that the
%   first_skip is required to be 0 in addition to the other conditions
%
dir_kind = cat(2, dir.kind);
dir_type = cat(2, dir.type);
dir_siz  = cat(2, dir.size);
dir_pos  = cat(2, dir.pos);
if all(dir_kind == FIFF.FIFF_DATA_BUFFER) && all(diff(dir_pos)==(dir_siz(1)+16)) && all(dir_siz==dir_siz(1)) && all(dir_type==dir_type(1)) && first_skip==0
  switch dir_type(1)
    case FIFF.FIFFT_DAU_PACK16
      nsamp = dir_siz(1)/(2*nchan);
    case FIFF.FIFFT_SHORT
      nsamp = dir_siz(1)/(2*nchan);
    case FIFF.FIFFT_FLOAT
      nsamp = dir_siz(1)/(4*nchan);
    case FIFF.FIFFT_INT
      nsamp = dir_siz(1)/(4*nchan);
    case FIFF.FIFFT_DOUBLE
      nsamp = dir_siz(1)/(8*nchan);
    case FIFF.FIFFT_COMPLEX_FLOAT
      nsamp = dir_siz(1)/(8*nchan);
    case FIFF.FIFFT_COMPLEX_DOUBLE
      nsamp = dir_siz(1)/(16*nchan);
    otherwise
      fclose(fid);
      error(me,'Cannot handle data buffers of type %d',dir_type(1));
  end
  nsamp = double(nsamp);
  
  firstsamp = cumsum([first_samp nsamp(ones(1,nent-first))]); 
  lastsamp  = firstsamp + nsamp - 1;
  
  % organize this slightly differently than in the original way. This is
  % much faster, and should not be too big of a problem, provided the
  % downstream reading function also knows how to interpret both
  % representations
  rawdir.first = firstsamp;
  rawdir.last  = lastsamp;
  rawdir.pos   = dir_pos;
  rawdir.type  = dir_type(1);
  rawdir.kind  = dir_kind(1);
  rawdir.size  = dir_siz(1);
  rawdir.nsamp = nsamp;
  
  data.first_samp = firstsamp(1);
  data.last_samp  = lastsamp(end);
  data.fastread   = true; % add flag that can be used later for fast reading
else
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
        case FIFF.FIFFT_DOUBLE
          nsamp = ent.size/(8*nchan);
        case FIFF.FIFFT_COMPLEX_FLOAT
          nsamp = ent.size/(8*nchan);
        case FIFF.FIFFT_COMPLEX_DOUBLE
          nsamp = ent.size/(16*nchan);
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
end

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

end
