function [data] = fiff_read_evoked_all(fname)
%
% [data] = fiff_read_evoked_all(fname)
%
% Read all evoked data set (averages only)
%


%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.6  2009/03/30 11:37:37  msh
%   Added copying of measurement info blocks from the original like in mne_browse_raw
%
%   Revision 1.5  2008/03/28 13:35:51  msh
%   Accommodate evoked files from xplotter, which may have channel information
%   local to the evoked block
%
%   Revision 1.4  2007/11/13 10:55:32  msh
%   Specify the second argument to all calls to the exist function
%
%   Revision 1.3  2006/10/14 11:47:32  msh
%   Omitted listing of nepoch from the diagnostic output
%
%   Revision 1.2  2006/10/14 11:39:58  msh
%   Fixed errors in fiff_read_evoked_all
%
%   Revision 1.1  2006/10/04 20:12:37  msh
%   Added fiff_read_evoked_all
%   Modified fiff_pick_channels_evoked to handle multiple data sets
%   Added bad channel writing to fiff_write_evoked
%
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_read_evoked_all';

if nargin ~= 1
    error(me,'Incorrect number of arguments');
end
%
%   Open the file
%
fprintf(1,'Reading %s ...\n',fname);
[ fid, tree ] = fiff_open(fname);
%
%   Read the measurement info
%
[ info, meas ] = fiff_read_meas_info(fid,tree);
info.filename = fname;
%
%   Locate the data of interest
%
processed = fiff_dir_tree_find(meas,FIFF.FIFFB_PROCESSED_DATA);
if length(processed) == 0
    fclose(fid);
    error(me,'Could not find processed data');
end
%
evoked = fiff_dir_tree_find(meas,FIFF.FIFFB_EVOKED);
if length(evoked) == 0
    fclose(fid);
    error(me,'Could not find evoked data');
end
%
%   Identify the aspects
%
naspect = 0;
is_smsh = [];
for k = 1:length(evoked)
    sets(k).aspects = fiff_dir_tree_find(evoked(k),FIFF.FIFFB_ASPECT);
    sets(k).naspect = length(sets(k).aspects);
    if sets(k).naspect > 0
        is_smsh = [ is_smsh zeros(1,sets(k).naspect) ];
        naspect = naspect + sets(k).naspect;
    end
    saspects  = fiff_dir_tree_find(evoked(k), FIFF.FIFFB_SMSH_ASPECT);
    nsaspects = length(saspects);
    if nsaspects > 0
        sets(k).naspect = sets(k).naspect + nsaspects;
        sets(k).aspects = [ sets(k).aspects saspects ];
        is_smsh = [ is_smsh ones(1,sets(k).naspect) ];
        naspect = naspect + nsaspects;
    end
end
if naspect == 0
    fclose(fid);
    error(me,'No evoked data aspects found');
end
fprintf(1,'\t%d evoked data sets containing a total of %d data aspects in %s\n',length(evoked),naspect,fname);
if ~isempty(info.comps)
    fprintf(1,'\t\t%d CTF compensation matrices available\n',length(info.comps));
end
%
%   Next locate the evoked data set
%
p = 1;
for e = 1:length(evoked)
    if e == 1
        %
        %   Find local channel information
        %
        my_evoked = evoked(e);
        nchan = 0;
        sfreq = -1;
        q = 0;
        for k = 1:my_evoked.nent
            kind = my_evoked.dir(k).kind;
            pos  = my_evoked.dir(k).pos;
            switch kind
                case FIFF.FIFF_NCHAN
                    tag = fiff_read_tag(fid,pos);
                    nchan = tag.data;
                case FIFF.FIFF_SFREQ
                    tag = fiff_read_tag(fid,pos);
                    sfreq = tag.data;
                case FIFF.FIFF_CH_INFO
                    q = q+1;
                    tag = fiff_read_tag(fid,pos);
                    chs(q) = tag.data;
            end
        end
        if nchan > 0
            if ~exist('chs','var')
                fclose(fid);
                error(me, ...
                    'Local channel information was not found when it was expected.');
            end
            if length(chs) ~= nchan
                fclose(fid);
                error(me, ...
                    'Number of channels and number of channel definitions are different');
            end
            info.chs   = chs;
            info.nchan = nchan;
            fprintf(1, ...
                '\tFound channel information in evoked data. nchan = %d\n',nchan);
            if sfreq > 0
                info.sfreq = sfreq;
            end
        end
    end
    for a = 1:sets(e).naspect
        my_evoked = evoked(e);
        my_aspect = sets(e).aspects(a);
        %
        %   Now find the data in the evoked block
        %
        for k = 1:my_evoked.nent
            kind = my_evoked.dir(k).kind;
            pos  = my_evoked.dir(k).pos;
            switch kind
                case FIFF.FIFF_COMMENT
                    tag = fiff_read_tag(fid,pos);
                    comment = tag.data;
                case FIFF.FIFF_FIRST_SAMPLE
                    tag = fiff_read_tag(fid,pos);
                    first = tag.data;
                case FIFF.FIFF_LAST_SAMPLE
                    tag = fiff_read_tag(fid,pos);
                    last = tag.data;
            end
        end
        if ~exist('comment','var')
            comment = 'No comment';
        end
        nsamp = last-first+1;
        %
        % Read the data in the aspect block
        %
        nepoch = 0;
        for k = 1:my_aspect.nent
            kind = my_aspect.dir(k).kind;
            pos  = my_aspect.dir(k).pos;
            switch kind
                case FIFF.FIFF_COMMENT
                    tag = fiff_read_tag(fid,pos);
                    comment = tag.data;
                case FIFF.FIFF_ASPECT_KIND
                    tag = fiff_read_tag(fid,pos);
                    aspect_kind = tag.data;
                case FIFF.FIFF_NAVE
                    tag = fiff_read_tag(fid,pos);
                    nave = tag.data;
                case FIFF.FIFF_EPOCH
                    nepoch = nepoch + 1;
                    tag = fiff_read_tag(fid,pos);
                    epoch(nepoch) = tag;
            end
        end
        if ~exist('nave','var')
            nave = 1;
        end
        if aspect_kind == FIFF.FIFFV_ASPECT_AVERAGE
            fprintf(1,'\tFound the data of interest:\n');
            fprintf(1,'\t\tt = %10.2f ... %10.2f ms (%s)\n',1000*first/info.sfreq,1000*last/info.sfreq,comment);
            fprintf(1,'\t\tnave = %d aspect type = %d\n',...
                nave,aspect_kind);
            if nepoch ~= 1 && nepoch ~= info.nchan
                fclose(fid);
                error(me,'Number of epoch tags is unreasonable (nepoch = %d nchan = %d)',nepoch,info.nchan);
            end
            %
            %       Put the old style epochs together
            %
            if nepoch == 1
                all = epoch(1).data;
            else
                all = epoch(1).data';
                for k = 2:nepoch
                    all = [ all ; epoch(k).data' ];
                end
            end
            if size(all,2) ~= nsamp
                fclose(fid);
                error(me,'Incorrect number of samples (%d instead of %d)',size(all,2),nsamp);
            end
            %
            %   Calibrate
            %
            for k = 1:info.nchan
                cals(k) = info.chs(k).cal;
            end
            all = diag(cals)*all;
            %
            %   Put it all together
            %
            one.aspect_kind = aspect_kind;
            one.is_smsh = is_smsh(p);
            if exist('nave','var')
                one.nave = nave;
            else
                one.nave  = 1;
            end
            one.first = first;
            one.last  = last;
            if exist('comment','var')
                one.comment = comment;
            end
            %
            %       Times for convenience and the actual epoch data
            %
            one.times = double(one.first:1:one.last)/info.sfreq;
            one.epochs = all;
            if p == 1
                data.info      = info;
                data.evoked(1) = one;
            else
                data.evoked(p) = one;
            end
            clear('one');
            p = p + 1;
        end
        clear('comment');
        clear('nave');
    end
end

fclose(fid);

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


