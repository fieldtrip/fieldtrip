function [leftmap,rightmap] = mne_read_morph_map(from,to,subjects_dir)
%
% [leftmap,rightmap] = mne_read_morph_map(from,to,subjects_dir)
%
% Read the morphing map from subject 'from' to subject 'to'.
% If subjects_dir is not specified, the SUBJECTS_DIR environment
% variable is used
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.1  2008/03/25 21:06:08  msh
%   Added mne_read_morph_map function
%
%
me='MNE:mne_read_morph_map';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin < 3
    subjects_dir=getenv('SUBJECTS_DIR');
    if isempty(subjects_dir)
        error(me,'SUBJECTS_DIR not set');
    end
end
if nargin < 2
    error(me,'Not enough input arguments');
end
%
%  Does the file exist
%
name = sprintf('%s/morph-maps/%s-%s-morph.fif',subjects_dir,from,to);
if ~exist(name,'file')
    name = sprintf('%s/morph-maps/%s-%s-morph.fif',subjects_dir,to,from);
    if ~exist(name,'file')
        error(me,'The requested morph map does not exist');
    end
end
%
%    Open it
%
[fid,tree] = fiff_open(name);
%
%   Locate all maps
%
maps = fiff_dir_tree_find(tree,FIFF.FIFFB_MNE_MORPH_MAP);
if isempty(maps)
    fclose(fid);
    error(me,'Morphing map data not found');
end
%
%   Find the correct ones
%
for k = 1:length(maps)
    tag = find_tag(maps(k),FIFF.FIFF_MNE_MORPH_MAP_FROM);
    if strcmp(tag.data,from)
        tag = find_tag(maps(k),FIFF.FIFF_MNE_MORPH_MAP_TO);
        if strcmp(tag.data,to)
            %
            %  Names match: which hemishere is this?
            %
            tag = find_tag(maps(k),FIFF.FIFF_MNE_HEMI);
            if tag.data == FIFF.FIFFV_MNE_SURF_LEFT_HEMI
                tag     = find_tag(maps(k),FIFF.FIFF_MNE_MORPH_MAP);
                leftmap = tag.data;
                fprintf(1,'\tLeft-hemisphere map read.\n');
            elseif tag.data == FIFF.FIFFV_MNE_SURF_RIGHT_HEMI
                tag      = find_tag(maps(k),FIFF.FIFF_MNE_MORPH_MAP);
                rightmap = tag.data;
                fprintf(1,'\tRight-hemisphere map read.\n');
            end
        end
    end
end
fclose(fid);
if ~exist('leftmap')
    error(me,'Left hemisphere map not found in %s',name);
end
if ~exist('rightmap')
    error(me,'Left hemisphere map not found in %s',name);
end

return;

    function [tag] = find_tag(node,findkind)
        
        for p = 1:node.nent
            if node.dir(p).kind == findkind
                tag = fiff_read_tag(fid,node.dir(p).pos);
                return;
            end
        end
        tag = [];
        return;
    end

end
