function [datout, S] = smudge(datin, tri, niter)

% SMUDGE(DATIN, TRI) computes a smudged version of the input data datain,
% given a triangulation tri. The algorithm is according to what is in
% MNE-Suite, documented in chapter 8.3

if nargin<3,
  niter = 1;
end

for k = 1:niter
  [tmp, Stmp] = do_smudge(datin, tri);
  if k==1,
    S      = Stmp;
  else
    S      = Stmp*S;
  end
  datout = tmp;
  datin  = tmp;
end

function [datout, S] = do_smudge(datin, tri)

% number of points
npnt = numel(datin);

% non-zero points
nz  = find(datin);

% triangles containing non-zero points
nzt = ismember(tri, nz);
tnz = sum(nzt, 2)>0 & sum(nzt,2)<3;

% points with non-zero neighbours
nzn = tri(tnz, :);
nzt = nzt(tnz, :);

vec1 = zeros(size(nzn,1)*2,1);
vec2 = zeros(size(nzn,1)*2,1);
for k = 1:size(nzn,1)
  indx0 = (k-1)*2+(1:2);
  indx1 = nzn(k, nzt(k,:));
  indx2 = nzn(k,~nzt(k,:));
  
  vec1(indx0) = indx1;
  vec2(indx0) = indx2;
  
end

vecx = sortrows([vec1 vec2]);
vecy = sortrows([vec2 vec1]);
clear vec1 vec2;

% in a closed surface the edges occur in doubles
vecx = vecx(1:2:end,:);
vecy = vecy(1:2:end,:);

tmp   = diff([0;vecx(:,1)])>0;
begix = find([tmp;1]==1); 

uvecy = unique(vecy(:,1));
tmp   = diff([0;vecy(:,1)])>0;
nix   = diff(find([tmp;1]==1));
nix(end+1:numel(uvecy)) = 1;

tmpval = zeros(npnt,1);
for k = 1:numel(uvecy)
  tmpval(uvecy(k)) = nix(k);
end

val = zeros(size(vecx,1),1);
for k = 1:numel(val)
  val(k) = 1./tmpval(vecx(k,2));
end

% allocate memory
S    = spalloc(npnt, npnt, numel(val)+numel(nz)); 
for k = 1:numel(begix)-1
  indx1 = vecx(begix(k),1);
  indx2 = begix(k):(begix(k+1)-1);
  indx3 = vecx(indx2, 2);
  S(indx3, indx1) = val(indx2);  
  %FIXME this can be speeded up by using the sparse command
end
S = S + spdiags(datin(:)>0, 0, npnt, npnt);

datout = S*datin(:);
