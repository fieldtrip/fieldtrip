function [data_sets] = fiff_find_evoked(fname)
%
% [data_sets] = fiff_find_evoked(fname)
%
% Find all evoked data sets in a fif file and create a list of descriptors
%


%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me = 'MNE:fiff_find_evoked';

%
%   Open the file
%
[ fid, tree ] = fiff_open(fname);
data_sets = struct('comment',{},'aspect_kind',{},'aspect_name',{});
%
%   Find all evoked data sets
%
evoked = fiff_dir_tree_find(tree, FIFF.FIFFB_EVOKED);
if length(evoked) == 0
    fclose(fid);
    return
end
%
%   Identify the aspects
%
naspect = 0;
for k = 1:length(evoked)
    sets(k).aspects  = fiff_dir_tree_find(evoked(k), FIFF.FIFFB_ASPECT);
    sets(k).naspect = length(sets(k).aspects);
    naspect = naspect + sets(k).naspect;
end
%
%   Collect the desired information
%
count = 1;
for k = 1:length(evoked)
    evoked_comment = find_tag(evoked(k), FIFF.FIFF_COMMENT);
    for a = 1:sets(k).naspect
        aspect_comment = find_tag(sets(k).aspects(a), FIFF.FIFF_COMMENT);
        aspect_kind    = find_tag(sets(k).aspects(a), FIFF.FIFF_ASPECT_KIND);
        if ~isempty(aspect_comment)
            data_sets(count).comment = aspect_comment.data;
        elseif ~isempty(evoked_comment)
            data_sets(count).comment = evoked_comment.data;
        else
            data_sets(count).comment = 'No comment';
        end
        if ~isempty(aspect_kind)
            data_sets(count).aspect_kind = aspect_kind.data;
        else
            data_sets(count).aspect_kind = -1;
        end
        switch data_sets(count).aspect_kind
        case FIFF.FIFFV_ASPECT_AVERAGE
            data_sets(count).aspect_name = 'Average';
        case FIFF.FIFFV_ASPECT_STD_ERR
            data_sets(count).aspect_name = 'Standard error';
        case FIFF.FIFFV_ASPECT_SINGLE
            data_sets(count).aspect_name = 'Single';
        case FIFF.FIFFV_ASPECT_SUBAVERAGE
            data_sets(count).aspect_name = 'Subaverage';
        case FIFF.FIFFV_ASPECT_ALTAVERAGE
            data_sets(count).aspect_name = 'Alt. subaverage';
        case FIFF.FIFFV_ASPECT_SAMPLE
            data_sets(count).aspect_name = 'Sample';
        case FIFF.FIFFV_ASPECT_POWER_DENSITY
            data_sets(count).aspect_name = 'Power density spectrum';
        case FIFF.FIFFV_ASPECT_DIPOLE_WAVE
            data_sets(count).aspect_name = 'Dipole source waveform';
            otherwise
            data_sets(count).aspect_name = 'Unknown';
        end
        count = count + 1;
    end
end

fclose(fid);

return

function [tag] = find_tag(node, findkind)

for p = 1:node.nent
    kind = node.dir(p).kind;
    pos  = node.dir(p).pos;
    if kind == findkind
        tag = fiff_read_tag(fid, pos);
        return
    end
end
tag = [];
return
end

end
