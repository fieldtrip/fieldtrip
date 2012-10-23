function [datout, S] = smudge(datin, tri, niter, threshold)

% SMUDGE(DATIN, TRI) computes a smudged version of the input data datain,
% given a triangulation tri. The algorithm is according to what is in
% MNE-Suite, documented in chapter 8.3

if nargin<3 || isempty(niter),
  niter = 1;
end

if nargin<4
  threshold = 0;
end

for k = 1:niter
  [tmp, Stmp] = do_smudge(datin, tri, threshold);
  if k==1,
    S      = Stmp;
  else
    S      = Stmp*S;
  end
  datout = tmp;
  datin  = tmp;
end

function [datout, S] = do_smudge(datin, tri, threshold)

% number of points
npnt = numel(datin);

% non-zero points
nz  = find(datin>threshold);

% triangles containing 1 or 2 non-zero points 
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

% matrices that define edges which has one of the members 0
% sorted according to the left or the right column respectively
vecx = sortrows([vec1 vec2]); % having the non-zero vertex upfront, sorted accordingly
vecy = sortrows([vec2 vec1]); % having the zero vertex upfront, sorted accordingly
clear vec1 vec2;

% in a closed surface the edges occur in doubles
vecx = vecx(1:2:end,:);
vecy = vecy(1:2:end,:);

% vecx now has in the first column the column indices for matrix S
% vecx now has in the second column the row indices for matrix S

% reconstruct the value that has to be put into the matrix
[uval,i1,i2] = unique(vecy(:,1));
tmp   = diff([0;vecy(:,1)])>0;
nix   = diff(find([tmp;1]==1));
nix(end+1:numel(uval)) = 1;

[~,i1,i2] = unique(vecx(:,2));
val  = 1./nix(i2);

S = sparse(vecx(:,2),vecx(:,1),val,npnt,npnt);
S = S + spdiags(datin(:)>threshold, 0, npnt, npnt);

datout = S*datin(:);
