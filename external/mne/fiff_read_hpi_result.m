function [ res ] = fiff_read_hpi_result(name)

%
% [ res ] = fiff_read_hpi_result(name)
%
% Read the HPI result block from a measurement file
%



global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_read_hpi_result';

if nargin ~= 1
    error(me,'Incorrect number of arguments');
end

[ fid, tree ] = fiff_open(name);

hpi = fiff_dir_tree_find(tree,FIFF.FIFFB_HPI_RESULT);

k = 0;
for p = 1:hpi.nent
    kind = hpi.dir(p).kind;
    pos  = hpi.dir(p).pos;
    switch kind
        case FIFF.FIFF_DIG_POINT
            k = k + 1;
            tag = fiff_read_tag(fid,pos);
            res.dig(k) = tag.data;
            res.dig(k).coord_frame = FIFF.FIFFV_COORD_DEVICE;
        case FIFF.FIFF_HPI_DIGITIZATION_ORDER
            tag = fiff_read_tag(fid,pos);
            res.dig_order = tag.data;
        case FIFF.FIFF_HPI_FIT_GOODNESS
            tag = fiff_read_tag(fid,pos);
            res.goodness = tag.data;
        case FIFF.FIFF_HPI_COILS_USED
            tag = fiff_read_tag(fid,pos);
            res.used = tag.data;
        case FIFF.FIFF_HPI_FIT_GOOD_LIMIT
            tag = fiff_read_tag(fid,pos);
            res.accept_limit = tag.data;
        case FIFF.FIFF_HPI_FIT_DIST_LIMIT
            tag = fiff_read_tag(fid,pos);
            res.dist_limit = tag.data;
        case FIFF.FIFF_HPI_FIT_ACCEPT
            tag = fiff_read_tag(fid,pos);
            res.accepted = tag.data;
        case FIFF.FIFF_HPI_COIL_MOMENTS
            tag = fiff_read_tag(fid,pos);
            res.moments = tag.data;
        case FIFF.FIFF_COORD_TRANS
            tag = fiff_read_tag(fid,pos);
            res.trans = tag.data;
        end
end
res.coord_frame = FIFF.FIFFV_COORD_DEVICE;

fclose(fid);
        
return;

end

