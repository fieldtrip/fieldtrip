function [track] = read_trk(filename, refine)

%read TrackVis .trk format data
% fillPath: filename of track to read.
%for format details http://www.trackvis.org/docs/?subsect=fileformat


fid = fopen(filename);
hdr = get_hdr(fid);
if hdr.hdr_size~=1000
  fclose(fid);
  fid = fopen(filePath, 'r', 'b'); % Big endian for old PPCs
  hdr = get_hdr(fid);
  if(hdr.hdr_size ~= 1000)
    error('Header length is not 1000, file may be corrupted');
  end
end

% go to the first data point and read the data as floats
fseek(fid, 1000, -1);
data = fread(fid, inf, 'float32');

% go to the first data point and read the data as integers
fseek(fid, 1000, -1); 
points = fread(fid, inf, 'int32');

fclose(fid);

dat  = cell(hdr.n_count,1);
npos = zeros(hdr.n_count,1);
offset = 0;
for k = 1:hdr.n_count
  npos(k) = points(offset+1);
  offset  = offset+1;
  dat{k}  = reshape(data(offset+(1:(3+hdr.n_scalars)*npos(k)),:),[],npos(k))';
  offset  = offset+(3+hdr.n_scalars)*npos(k)+hdr.n_properties;
end
pos = cat(1,dat{:});
pos = pos*hdr.vox_to_ras(1:3,1:3)' + hdr.vox_to_ras(1:3,4*ones(1,size(pos,1)))';

% convert the tracks to lines
maxnpos = max(npos);
cumnpos = cumsum([0;npos]);
lines   = zeros(maxnpos*2, numel(npos))+nan;
for k = 1:numel(npos)
  lines(1:npos(k),k) = cumnpos(k)+(1:npos(k));
  lines(npos(k)+(1:npos(k)),k) = (cumnpos(k)+npos(k)):-1:(cumnpos(k)+1);
end
track.line = lines';
track.pos  = pos;
track.hdr  = hdr;


function hdr = get_hdr(fid)

hdr.id_string                 = fread(fid, 6, '*char')';
hdr.dim                       = fread(fid, 3, 'short')';
hdr.voxel_size                = fread(fid, 3, 'float')';
hdr.origin                    = fread(fid, 3, 'float')';
hdr.n_scalars                 = fread(fid, 1, 'short')';
hdr.scalar_name               = fread(fid, [20,10], '*char')';
hdr.n_properties              = fread(fid, 1, 'short')';
hdr.property_name             = fread(fid, [20,10], '*char')';
hdr.vox_to_ras                = fread(fid, [4,4], 'float')';
hdr.reserved                  = fread(fid, 444, '*char');
hdr.voxel_order               = fread(fid, 4, '*char')';
hdr.pad2                      = fread(fid, 4, '*char')';
hdr.image_orientation_patient = fread(fid, 6, 'float')';
hdr.pad1                      = fread(fid, 2, '*char')';
hdr.invert_x                  = fread(fid, 1, 'uchar');
hdr.invert_y                  = fread(fid, 1, 'uchar');
hdr.invert_z                  = fread(fid, 1, 'uchar');
hdr.swap_xy                   = fread(fid, 1, 'uchar');
hdr.swap_yz                   = fread(fid, 1, 'uchar');
hdr.swap_zx                   = fread(fid, 1, 'uchar');
hdr.n_count                   = fread(fid, 1, 'int')';
hdr.version                   = fread(fid, 1, 'int')';
hdr.hdr_size                  = fread(fid, 1, 'int')';
%end getHeader()