function data = read_tck(filename)

% READ_TCK reads tractography information from an mrtrix-generated .tck
% file. Requires the matlab functions from mrtrix.

tmp  = read_mrtrix_tracks(filename);
npos = cellfun(@numel, tmp.data)./3;
maxnpos = max(npos);
lines   = zeros(2*maxnpos, numel(npos))+nan;
offset  = 0;
for k = 1:numel(npos)
  indx1 = 1:npos(k)*2;
  indx2 = [1:npos(k) npos(k):-1:1]+offset;
  lines(indx1,k) = indx2;
  offset = max(indx2);
end

data.pos  = cat(1, tmp.data{:});
data.line = lines.';
data.hdr  = rmfield(tmp, 'data');
